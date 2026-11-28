# Run one TABLES case.
#
#     julia --project run_case.jl flat_u00
#
# The case is one of flat_u00, flat_u10, ridge100_u00, ridge100_u10, ridge1000_u00, ridge1000_u10.
# The defaults follow the protocol at full size (512 x 1024 columns at 20 m, time step 0.2 s, 3 h).
# Environment variables override them, mainly to run reduced configurations for testing:
#
#     TABLES_ARCH               cpu or gpu (default: gpu)
#     TABLES_FLOAT_TYPE         Float64 or Float32 (default: Float64)
#     TABLES_NX, TABLES_NY      number of columns in x and y (default: 512, 1024)
#     TABLES_DX                 horizontal spacing in m (default: 20)
#     TABLES_DZ                 vertical spacing in m below TABLES_UNIFORM_HEIGHT (default: 10)
#     TABLES_UNIFORM_HEIGHT     height in m below which the vertical spacing is uniform (default: 3000)
#     TABLES_MAX_DZ             largest vertical spacing in m (default: 40)
#     TABLES_DT                 time step in s (default: 0.2)
#     TABLES_STOP_TIME          duration in s (default: 10800)
#     TABLES_SECTION_INTERVAL   time in s between outputs of y-averaged sections (default: 60)
#     TABLES_SLICE_INTERVAL     time in s between outputs of instantaneous x-z and x-y slices (default: 60)
#     TABLES_PROGRESS_INTERVAL  iterations between progress messages (default: 500)
#     TABLES_STATISTICS_WINDOW  length in s of the averaging window at the end of the run (default: 1800)
#     TABLES_STATISTICS_STRIDE  iterations between samples of the statistics (default: 1)
#     TABLES_OUTPUT_DIR         output directory (default: output)
#     TABLES_PICKUP             true to restart from the latest checkpoint (default: false)
#     TABLES_WALL_TIME_LIMIT    wall time in s after which the run stops and writes a checkpoint, to be
#                               continued with TABLES_PICKUP=true (default: no limit)
#
# jobs/submit_case.sh runs a case on a Casper GPU as a chain of such runs.

include(joinpath(@__DIR__, "tables.jl"))
include(joinpath(@__DIR__, "protocol_output.jl"))

setting(T, name, default) = haskey(ENV, name) ? parse(T, ENV[name]) : default

case_name = Symbol(ARGS[1])

arch_name = lowercase(get(ENV, "TABLES_ARCH", "gpu"))
if arch_name == "gpu"
    using CUDA
    arch = GPU()
else
    arch = CPU()
end

FT = get(ENV, "TABLES_FLOAT_TYPE", "Float64") == "Float32" ? Float32 : Float64
Oceananigans.defaults.FloatType = FT

Nx = setting(Int, "TABLES_NX", 512)
Ny = setting(Int, "TABLES_NY", 1024)
Δx = setting(Float64, "TABLES_DX", 20)
Δt = setting(Float64, "TABLES_DT", 0.2)
stop_time = setting(Float64, "TABLES_STOP_TIME", 3hours)
output_dir = get(ENV, "TABLES_OUTPUT_DIR", "output")

z_faces = vertical_faces(spacing = setting(Float64, "TABLES_DZ", 10),
                         uniform_height = setting(Float64, "TABLES_UNIFORM_HEIGHT", 3000),
                         maximum_spacing = setting(Float64, "TABLES_MAX_DZ", 40))

case = tables_case(case_name, FT)
grid = tables_grid(arch, FT, case; Nx, Ny, Δx, z_faces)

# In single precision the increments of a 0.2 s step are close to the resolution of ρθ.
# Stochastic rounding removes the resulting bias (see "Round-off in reduced precision" in the Breeze docs).
timestepper_kwargs = FT === Float32 ? (; rounding = StochasticRounding()) : NamedTuple()

model = tables_model(grid, case; timestepper_kwargs)
set_tables_initial_condition!(model, case)

@info "TABLES case $case_name" grid model
@info @sprintf("Pressure imbalance of the initial state: %.2e Pa", maximum_pressure_imbalance(model))

simulation = tables_simulation(model; Δt, stop_time,
                               progress_interval = setting(Int, "TABLES_PROGRESS_INTERVAL", 500))

# A run with a wall time limit stops by itself before its job is killed
simulation.wall_time_limit = setting(Float64, "TABLES_WALL_TIME_LIMIT", Inf)

#####
##### Monitoring output, for checking that a case runs correctly
#####

pickup = get(ENV, "TABLES_PICKUP", "false") == "true"

# `run!` only warns when it finds no checkpoint to pick up, and then starts from the initial condition
if pickup && !(isdir(output_dir) && any(startswith("$(case_name)_checkpoint_iteration"), readdir(output_dir)))
    error("TABLES_PICKUP=true, but there is no checkpoint of $case_name in $output_dir")
end

iterations(interval) = max(1, round(Int, interval / Δt))

u, v, w = model.velocities
θ = model.formulation.potential_temperature

outputs = (; u, v, w, θ, ww = w^2, wθ = w * θ)

# Averages along the ridge, on model levels
sections = NamedTuple(name => Average(@at((Center, Center, Center), outputs[name]), dims = 2)
                      for name in keys(outputs))

section_interval = setting(Float64, "TABLES_SECTION_INTERVAL", 1minute)

simulation.output_writers[:sections] = JLD2Writer(model, sections;
                                                  dir = output_dir,
                                                  filename = "$(case_name)_sections",
                                                  schedule = IterationInterval(iterations(section_interval)),
                                                  overwrite_files = !pickup)

# Instantaneous slices on model levels (see animate.jl): the x-z slice in the middle of the domain
# and the x-y slice on the level whose w is closest to 500 m. There u, v and θ are half a cell higher.
slices = (; u, v, w, θ)
slice_interval = setting(Float64, "TABLES_SLICE_INTERVAL", 1minute)

simulation.output_writers[:slices] = JLD2Writer(model, slices;
                                                dir = output_dir,
                                                filename = "$(case_name)_slices",
                                                indices = (:, Ny ÷ 2, :),
                                                schedule = IterationInterval(iterations(slice_interval)),
                                                overwrite_files = !pickup)

simulation.output_writers[:xy_slices] = JLD2Writer(model, slices;
                                                   dir = output_dir,
                                                   filename = "$(case_name)_xy_slices",
                                                   indices = (:, :, argmin(abs.(z_faces .- 500))),
                                                   schedule = IterationInterval(iterations(slice_interval)),
                                                   overwrite_files = !pickup)

# Checkpoints every 30 minutes. The one at 2.5 h is the start of the statistics window.
simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                        dir = output_dir,
                                                        prefix = "$(case_name)_checkpoint",
                                                        schedule = IterationInterval(iterations(30minutes)))

#####
##### Output requested by the protocol (see protocol_output.jl)
#####

# Statistics averaged over y and over the last 30 minutes (Sections 3 and 4). The callback has a
# name so that a restart finds its sums in the checkpoint.
statistics_window = setting(Float64, "TABLES_STATISTICS_WINDOW", 30minutes)
statistics = TurbulenceStatistics(grid; window_start = stop_time - statistics_window)
add_callback!(simulation, statistics, IterationInterval(setting(Int, "TABLES_STATISTICS_STRIDE", 1)); name = :statistics)

# A run that the wall time limit can stop also writes a checkpoint where it stops
run!(simulation; pickup, checkpoint_at_end = isfinite(simulation.wall_time_limit))

if time(simulation) < simulation.stop_time
    @info "Stopped by the wall time limit at t = $(prettytime(simulation)). Continue with TABLES_PICKUP=true."
else
    statistics_file = write_statistics(joinpath(output_dir, "$(case_name)_statistics.nc"), statistics, model; case_name)
    @info @sprintf("Wrote %s with %d samples from t = %s to t = %s", statistics_file, statistics.samples[],
                   prettytime(statistics.first_sample_time[]), prettytime(statistics.last_sample_time[]))

    # Snapshot of the final state (Section 3.1)
    snapshot_file = write_snapshot(joinpath(output_dir, "$(case_name)_snapshot.nc"), model; case_name)
    @info "Wrote $snapshot_file"

    # This file tells jobs/run_case.pbs that the case is complete
    touch(joinpath(output_dir, "$(case_name).done"))
end
