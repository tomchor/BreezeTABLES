#!/bin/bash
# Submit a TABLES case to a Casper GPU as a chain of jobs, followed by the job that animates it.
#
#     jobs/submit_case.sh flat_u00            # the protocol at full size
#     jobs/submit_case.sh flat_u00 RUN=test WALLTIME=00:30:00 TABLES_WALL_TIME_LIMIT=600 \
#                         TABLES_NX=64 TABLES_NY=128 TABLES_STOP_TIME=600 TABLES_STATISTICS_WINDOW=300
#
# Each job (run_case.pbs) is one segment of the run: it stops before its wall time, writes a
# checkpoint and submits the next segment, which continues from it. The last one writes the
# statistics and the snapshot and submits the animation (animate.pbs). Everything goes to the run
# directory /glade/derecho/scratch/$USER/BreezeTABLES/$RUN.
#
# The jobs use a copy of the scripts and of the Julia environment that is taken here, in code/ of
# the run directory, so that later changes to the repository do not reach a run that is under way.
#
# Settings, given as NAME=value after the case:
#   RUN            name of the run directory (default: the case)
#   GPU_TYPE       cc80 (default) for an A100 or an H100, the GPUs with the 80 GB that the full size
#                  needs in double precision; or a100, h100, v100, ...
#   WALLTIME       wall time of a segment (default: 12:00:00)
#   MAX_SEGMENTS   the chain stops after this many segments, complete or not (default: 8)
#   TABLES_*       settings of run_case.jl. TABLES_WALL_TIME_LIMIT defaults to WALLTIME minus one hour,
#                  which leaves time to start Julia and to write the last checkpoint and the final files.
#
# The settings are written to settings.env in the run directory and apply to every segment.
#
# To continue a run by hand, for example after a job failed, and to animate a run by hand,
# complete or not:
#
#     jobs/submit_case.sh --continue <run directory> <segment>
#     jobs/submit_case.sh --animate <run directory>
set -euo pipefail

submit_segment() { # <run directory> <segment>
    source "$1/settings.env"
    qsub -N "$RUN-$2" \
         -l "select=1:ncpus=8:ngpus=1:gpu_type=$GPU_TYPE:mem=64GB" -l "walltime=$WALLTIME" \
         -o "$1/logs/segment$2.pbs.log" \
         -v "RUN_DIR=$1,SEGMENT=$2" \
         "$CODE/jobs/run_case.pbs"
}

submit_animation() { # <run directory>
    source "$1/settings.env"
    qsub -N "$RUN-movie" -o "$1/logs/animation.pbs.log" -v "RUN_DIR=$1" "$CODE/jobs/animate.pbs"
}

case ${1:-} in
    --continue) submit_segment "$2" "$3"; exit ;;
    --animate) submit_animation "$2"; exit ;;
esac

CASE=$1
shift

RUN=$CASE
GPU_TYPE=cc80
WALLTIME=12:00:00
MAX_SEGMENTS=8
tables_settings=()

for setting in "$@"; do
    case $setting in
        RUN=* | GPU_TYPE=* | WALLTIME=* | MAX_SEGMENTS=*) declare "$setting" ;;
        TABLES_*=*) tables_settings+=("$setting") ;;
        *) echo "Unknown setting: $setting"; exit 1 ;;
    esac
done

repository=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
run_directory=/glade/derecho/scratch/$USER/BreezeTABLES/$RUN

if [ -e "$run_directory/settings.env" ]; then
    echo "$run_directory already holds a run"
    exit 1
fi

if [[ " ${tables_settings[*]} " != *" TABLES_WALL_TIME_LIMIT="* ]]; then
    IFS=: read -r hours minutes seconds <<< "$WALLTIME"
    wall_time_limit=$((10#$hours * 3600 + 10#$minutes * 60 + 10#$seconds - 3600))

    if [ "$wall_time_limit" -le 0 ]; then
        echo "WALLTIME=$WALLTIME leaves no time for the default wall time limit: set TABLES_WALL_TIME_LIMIT"
        exit 1
    fi

    tables_settings+=("TABLES_WALL_TIME_LIMIT=$wall_time_limit")
fi

mkdir -p "$run_directory/logs" "$run_directory/code"

cp "$repository"/*.jl "$repository"/*.toml "$run_directory/code"
cp -r "$repository/jobs" "$run_directory/code"
echo "$(git -C "$repository" describe --always --dirty) of $repository" > "$run_directory/code/source"

{
    echo "RUN=$RUN"
    echo "CASE=$CASE"
    echo "CODE=$run_directory/code"
    echo "GPU_TYPE=$GPU_TYPE"
    echo "WALLTIME=$WALLTIME"
    echo "MAX_SEGMENTS=$MAX_SEGMENTS"
    for setting in "${tables_settings[@]}"; do
        echo "$setting"
    done
} > "$run_directory/settings.env"

submit_segment "$run_directory" 1
