# Output requested by the TABLES protocol (Sections 3 to 5 of Rauchoecker_QJRMS_TABLES.pdf).
#
# This file has the parts that must run with the simulation:
#
#   - the turbulence statistics, accumulated during the last 30 minutes (Sections 3 and 4), and
#   - the snapshot of the final state (Section 3.1).
#
# Both are written on model levels. Interpolating them to the common grid of Section 5, computing
# the spectra from the snapshot, and writing the files under the names the protocol asks for are
# post-processing steps that are not implemented yet.
#
# Not available without an explicit closure: the subfilter-scale stresses and heat fluxes
# (`*_sfs`) and the dissipation rates (`eps`, `eps_t`).

using Breeze
using Oceananigans
using Oceananigans: prognostic_state, restore_prognostic_state!
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: znode
using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑzᵃᵃᶜ, Δzᶜᶜᶜ
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index
using NCDatasets

#####
##### Statistics (Sections 3 and 4)
#####
##### The protocol allows computing the statistics online from the relation
#####
#####   ⟨a′b′⟩ = ⟨ab⟩ - ⟨a⟩⟨b⟩,
#####
##### where ⟨⟩ is the average over time and over y. We accumulate the sums of the raw moments
##### listed below at every cell center and form the central moments at the end. The sums are
##### kept in double precision whatever the precision of the model, because the central moments
##### are small differences of large numbers (θ′² is 10⁻⁷ of θ²).
#####

const moment_names = (:u, :v, :w, :θ, :p, :ρ,
                      :uu, :uv, :uw, :vv, :vw, :ww,
                      :uθ, :vθ, :wθ, :θθ,
                      :up, :vp, :wp,
                      :uuu, :uuv, :uuw, :vvu, :vvv, :vvw, :wwu, :wwv, :www,
                      :uuθ, :vvθ, :wwθ)

const number_of_moments = length(moment_names)

@inline moments(u, v, w, θ, p, ρ) = (u, v, w, θ, p, ρ,
                                      u * u, u * v, u * w, v * v, v * w, w * w,
                                      u * θ, v * θ, w * θ, θ * θ,
                                      u * p, v * p, w * p,
                                      u * u * u, u * u * v, u * u * w,
                                      v * v * u, v * v * v, v * v * w,
                                      w * w * u, w * w * v, w * w * w,
                                      u * u * θ, v * v * θ, w * w * θ)

@inline add_tuples(a::NTuple{N}, b::NTuple{N}) where N = ntuple(n -> @inbounds(a[n] + b[n]), Val(N))

# One work item per (i, k). Each sums the moments along y and adds the result to the running sums.
@kernel function _accumulate_moments!(moment_sums, friction_velocity_sum, grid, u, v, w, θ, p, ρ, Ny, wall)
    i, k = @index(Global, NTuple)

    sums = ntuple(n -> 0.0, Val(number_of_moments))
    friction_velocity = 0.0

    @inbounds for j in 1:Ny
        uᶜ = Float64(ℑxᶜᵃᵃ(i, j, k, grid, u))
        vᶜ = Float64(ℑyᵃᶜᵃ(i, j, k, grid, v))
        wᶜ = Float64(ℑzᵃᵃᶜ(i, j, k, grid, w))

        sums = add_tuples(sums, moments(uᶜ, vᶜ, wᶜ, Float64(θ[i, j, k]), Float64(p[i, j, k]), Float64(ρ[i, j, k])))

        # Friction velocity as defined in Section 3.1, from the wind at the lowest level
        if k == 1
            z₁ = Float64(Δzᶜᶜᶜ(i, j, 1, grid)) / 2
            friction_velocity += wall.von_karman_constant * sqrt(uᶜ^2 + vᶜ^2) / log(z₁ / wall.roughness_length)
        end
    end

    @inbounds for n in 1:number_of_moments
        moment_sums[i, k, n] += sums[n]
    end

    if k == 1
        @inbounds friction_velocity_sum[i] += friction_velocity
    end
end

"""
    TurbulenceStatistics(grid; window_start, roughness_length = 0.1, von_karman_constant = 0.4)

Accumulator of the turbulence statistics of the protocol, averaged over y and over the simulation
time after `window_start` (in seconds). Add it to a simulation as a callback:

    statistics = TurbulenceStatistics(grid; window_start = 2.5hours)
    add_callback!(simulation, statistics, IterationInterval(1))

Every call after `window_start` takes one sample of the model state. All fields are sampled at
cell centers, with the velocity components interpolated there. `write_statistics` forms the
statistics and writes them.

The accumulated sums are saved in the checkpoints of the simulation and restored from them, so a
simulation can be interrupted and restarted inside the averaging window.
"""
struct TurbulenceStatistics{A, V, W}
    moment_sums :: A            # (Nx, Nz, number_of_moments): sums over y and over the samples
    friction_velocity_sum :: V  # (Nx,)
    wall :: W
    window_start :: Float64
    samples :: Base.RefValue{Int}
    first_sample_time :: Base.RefValue{Float64}
    last_sample_time :: Base.RefValue{Float64}
end

function TurbulenceStatistics(grid; window_start, roughness_length = 0.1, von_karman_constant = 0.4)
    arch = architecture(grid)
    Nx, Ny, Nz = size(grid)
    moment_sums = on_architecture(arch, zeros(Float64, Nx, Nz, number_of_moments))
    friction_velocity_sum = on_architecture(arch, zeros(Float64, Nx))
    wall = (; von_karman_constant = Float64(von_karman_constant), roughness_length = Float64(roughness_length))

    return TurbulenceStatistics(moment_sums, friction_velocity_sum, wall, Float64(window_start),
                                Ref(0), Ref(NaN), Ref(NaN))
end

function (statistics::TurbulenceStatistics)(simulation)
    # Half a time step of tolerance, because the model time is a sum of floating-point steps
    time(simulation) ≥ statistics.window_start - simulation.Δt / 2 || return nothing

    # When rounding leaves a sliver of time before the stop time, the simulation skips that step
    # and calls back a second time at its last iteration. That state has already been sampled.
    time(simulation) < statistics.last_sample_time[] + simulation.Δt / 2 && return nothing

    sample!(statistics, simulation.model)
    return nothing
end

# What a checkpoint holds of the statistics
restart_state(statistics::TurbulenceStatistics) =
    (; statistics.moment_sums, statistics.friction_velocity_sum,
       statistics.samples, statistics.first_sample_time, statistics.last_sample_time)

Oceananigans.prognostic_state(statistics::TurbulenceStatistics) = prognostic_state(restart_state(statistics))
Oceananigans.restore_prognostic_state!(statistics::TurbulenceStatistics, from) = restore_prognostic_state!(restart_state(statistics), from)
Oceananigans.restore_prognostic_state!(::TurbulenceStatistics, ::Nothing) = nothing

"Add one sample of the state of `model` to `statistics`."
function sample!(statistics::TurbulenceStatistics, model)
    grid = model.grid
    arch = architecture(grid)
    Nx, Ny, Nz = size(grid)

    u, v, w = model.velocities
    θ = model.formulation.potential_temperature
    p = model.dynamics.pressure
    ρ = model.dynamics.total_density

    launch!(arch, grid, (Nx, Nz), _accumulate_moments!,
            statistics.moment_sums, statistics.friction_velocity_sum, grid, u, v, w, θ, p, ρ, Ny, statistics.wall)

    statistics.samples[] += 1
    isnan(statistics.first_sample_time[]) && (statistics.first_sample_time[] = model.clock.time)
    statistics.last_sample_time[] = model.clock.time

    return nothing
end

covariance(ab, a, b) = ab .- a .* b

# ⟨a′a′b′⟩ from the raw moments ⟨aab⟩, ⟨aa⟩, ⟨ab⟩, ⟨a⟩ and ⟨b⟩
third_moment(aab, aa, ab, a, b) = aab .- 2 .* a .* ab .- b .* aa .+ 2 .* a.^2 .* b

"""
    protocol_statistics(statistics, grid)

Return a `NamedTuple` of `(Nx, Nz)` arrays with the statistics of Section 3 that have been
accumulated in `statistics`, under the variable names of the protocol, and the friction velocity
`ust` as a vector of length `Nx`.
"""
function protocol_statistics(statistics::TurbulenceStatistics, grid)
    samples = statistics.samples[]
    samples > 0 || error("No samples have been taken. Did the simulation reach the averaging window?")

    count = samples * size(grid, 2)
    means = Array(statistics.moment_sums) ./ count
    m = NamedTuple{moment_names}(Tuple(means[:, :, n] for n in 1:number_of_moments))

    fields = (; u_mean = m.u, v_mean = m.v, w_mean = m.w, t_mean = m.θ, p_mean = m.p, rho = m.ρ,

                upup_res = covariance(m.uu, m.u, m.u),
                upvp_res = covariance(m.uv, m.u, m.v),
                upwp_res = covariance(m.uw, m.u, m.w),
                vpvp_res = covariance(m.vv, m.v, m.v),
                vpwp_res = covariance(m.vw, m.v, m.w),
                wpwp_res = covariance(m.ww, m.w, m.w),

                tptp_res = covariance(m.θθ, m.θ, m.θ),
                uptp_res = covariance(m.uθ, m.u, m.θ),
                vptp_res = covariance(m.vθ, m.v, m.θ),
                wptp_res = covariance(m.wθ, m.w, m.θ),

                uppp = covariance(m.up, m.u, m.p),
                vppp = covariance(m.vp, m.v, m.p),
                wppp = covariance(m.wp, m.w, m.p),

                upupup = third_moment(m.uuu, m.uu, m.uu, m.u, m.u),
                upupvp = third_moment(m.uuv, m.uu, m.uv, m.u, m.v),
                upupwp = third_moment(m.uuw, m.uu, m.uw, m.u, m.w),
                vpvpup = third_moment(m.vvu, m.vv, m.uv, m.v, m.u),
                vpvpvp = third_moment(m.vvv, m.vv, m.vv, m.v, m.v),
                vpvpwp = third_moment(m.vvw, m.vv, m.vw, m.v, m.w),
                wpwpup = third_moment(m.wwu, m.ww, m.uw, m.w, m.u),
                wpwpvp = third_moment(m.wwv, m.ww, m.vw, m.w, m.v),
                wpwpwp = third_moment(m.www, m.ww, m.ww, m.w, m.w),

                upuptp = third_moment(m.uuθ, m.uu, m.uθ, m.u, m.θ),
                vpvptp = third_moment(m.vvθ, m.vv, m.vθ, m.v, m.θ),
                wpwptp = third_moment(m.wwθ, m.ww, m.wθ, m.w, m.θ))

    ust = Array(statistics.friction_velocity_sum) ./ count

    return fields, ust
end

const statistics_attributes = Dict(
    "u_mean" => ("mean wind component along x", "m s-1"),
    "v_mean" => ("mean wind component along y", "m s-1"),
    "w_mean" => ("mean vertical wind component", "m s-1"),
    "t_mean" => ("mean potential temperature", "K"),
    "p_mean" => ("mean pressure", "Pa"),
    "rho" => ("mean density", "kg m-3"),
    "upup_res" => ("resolved variance of u", "m2 s-2"),
    "upvp_res" => ("resolved covariance of u and v", "m2 s-2"),
    "upwp_res" => ("resolved covariance of u and w", "m2 s-2"),
    "vpvp_res" => ("resolved variance of v", "m2 s-2"),
    "vpwp_res" => ("resolved covariance of v and w", "m2 s-2"),
    "wpwp_res" => ("resolved variance of w", "m2 s-2"),
    "tptp_res" => ("resolved variance of potential temperature", "K2"),
    "uptp_res" => ("resolved potential temperature flux along x", "K m s-1"),
    "vptp_res" => ("resolved potential temperature flux along y", "K m s-1"),
    "wptp_res" => ("resolved vertical potential temperature flux", "K m s-1"),
    "uppp" => ("resolved covariance of u and pressure", "Pa m s-1"),
    "vppp" => ("resolved covariance of v and pressure", "Pa m s-1"),
    "wppp" => ("resolved covariance of w and pressure", "Pa m s-1"),
    "upupup" => ("resolved third moment u'u'u'", "m3 s-3"),
    "upupvp" => ("resolved third moment u'u'v'", "m3 s-3"),
    "upupwp" => ("resolved third moment u'u'w'", "m3 s-3"),
    "vpvpup" => ("resolved third moment v'v'u'", "m3 s-3"),
    "vpvpvp" => ("resolved third moment v'v'v'", "m3 s-3"),
    "vpvpwp" => ("resolved third moment v'v'w'", "m3 s-3"),
    "wpwpup" => ("resolved third moment w'w'u'", "m3 s-3"),
    "wpwpvp" => ("resolved third moment w'w'v'", "m3 s-3"),
    "wpwpwp" => ("resolved third moment w'w'w'", "m3 s-3"),
    "upuptp" => ("resolved third moment u'u'theta'", "K m2 s-2"),
    "vpvptp" => ("resolved third moment v'v'theta'", "K m2 s-2"),
    "wpwptp" => ("resolved third moment w'w'theta'", "K m2 s-2"))

#####
##### Coordinates of the model levels
#####

# The terrain does not vary along y, so heights are evaluated at j = 1.
function level_heights(grid, ℓx, ℓz)
    cpu_grid = on_architecture(CPU(), grid)
    Nx = size(cpu_grid, 1)
    Nz = ℓz isa Face ? size(cpu_grid, 3) + 1 : size(cpu_grid, 3)
    return [znode(i, 1, k, cpu_grid, ℓx, Center(), ℓz) for i in 1:Nx, k in 1:Nz]
end

surface_heights(grid, ℓx) = level_heights(grid, ℓx, Face())[:, 1]

function common_attributes(model, case_name)
    return Dict("title" => "TABLES case $case_name",
                "model" => "Breeze.jl $(pkgversion(Breeze)) with Oceananigans.jl $(pkgversion(Oceananigans))",
                "vertical_grid" => "model levels; `height` gives the height of every level above z = 0",
                "time" => model.clock.time)
end

#####
##### Files
#####

"""
    write_statistics(path, statistics, model; case_name)

Write the statistics accumulated in `statistics` to the NetCDF file `path`, with the variable names
of the protocol and dimensions `(x, z)`. The fields are on the model levels, at cell centers.
Following Section 5, `x` is the distance from the first cell center.
"""
function write_statistics(path, statistics::TurbulenceStatistics, model; case_name)
    grid = model.grid
    fields, ust = protocol_statistics(statistics, grid)

    x = Array(xnodes(grid, Center()))
    heights = level_heights(grid, Center(), Center())

    mkpath(dirname(path))
    isfile(path) && rm(path)

    NCDataset(path, "c") do dataset
        attributes = common_attributes(model, case_name)
        attributes["averaging"] = "over y and over time, from the raw moments accumulated at every sample"
        attributes["averaging_start"] = statistics.first_sample_time[]
        attributes["averaging_end"] = statistics.last_sample_time[]
        attributes["samples_in_time"] = statistics.samples[]
        attributes["missing_variables"] = "subfilter-scale fluxes (*_sfs) and dissipation rates (eps, eps_t): the model has no explicit closure"
        for (name, value) in attributes
            dataset.attrib[name] = value
        end

        defDim(dataset, "x", size(heights, 1))
        defDim(dataset, "z", size(heights, 2))

        defVar(dataset, "x", x .- minimum(x), ("x",); attrib = Dict("long_name" => "distance from the first cell center", "units" => "m"))
        defVar(dataset, "z", collect(1:size(heights, 2)), ("z",); attrib = Dict("long_name" => "model level"))
        defVar(dataset, "height", heights, ("x", "z"); attrib = Dict("long_name" => "height of the cell centers", "units" => "m"))
        defVar(dataset, "surface_height", surface_heights(grid, Center()), ("x",); attrib = Dict("long_name" => "terrain height", "units" => "m"))

        for (name, values) in pairs(fields)
            long_name, units = statistics_attributes[string(name)]
            defVar(dataset, string(name), values, ("x", "z"); attrib = Dict("long_name" => long_name, "units" => units))
        end

        defVar(dataset, "ust", ust, ("x",); attrib = Dict("long_name" => "friction velocity from the neutral log law at the lowest level", "units" => "m s-1"))
    end

    return path
end

"""
    write_snapshot(path, model; case_name)

Write the instantaneous `u`, `v`, `w`, potential temperature `t` and pressure `p` of `model` to the
NetCDF file `path` in single precision. The fields are on the model levels and at their own
locations on the staggered grid: `u` at `xf`, `v` at `yf`, `w` at `zf`, and `t` and `p` at cell centers.
"""
function write_snapshot(path, model; case_name)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    u, v, w = model.velocities
    θ = model.formulation.potential_temperature
    p = model.dynamics.pressure

    coordinate(nodes) = Array(nodes)[1:end]

    mkpath(dirname(path))
    isfile(path) && rm(path)

    NCDataset(path, "c") do dataset
        for (name, value) in common_attributes(model, case_name)
            dataset.attrib[name] = value
        end

        defDim(dataset, "x", Nx);  defDim(dataset, "xf", Nx)
        defDim(dataset, "y", Ny);  defDim(dataset, "yf", Ny)
        defDim(dataset, "z", Nz);  defDim(dataset, "zf", Nz + 1)

        meters = Dict("units" => "m")
        defVar(dataset, "x",  coordinate(xnodes(grid, Center())), ("x",);  attrib = meters)
        defVar(dataset, "xf", coordinate(xnodes(grid, Face())),   ("xf",); attrib = meters)
        defVar(dataset, "y",  coordinate(ynodes(grid, Center())), ("y",);  attrib = meters)
        defVar(dataset, "yf", coordinate(ynodes(grid, Face())),   ("yf",); attrib = meters)
        defVar(dataset, "z",  collect(1:Nz),     ("z",);  attrib = Dict("long_name" => "model level"))
        defVar(dataset, "zf", collect(1:Nz + 1), ("zf",); attrib = Dict("long_name" => "model level interface"))

        defVar(dataset, "height",   level_heights(grid, Center(), Center()), ("x", "z");  attrib = Dict("long_name" => "height of t and p", "units" => "m"))
        defVar(dataset, "height_u", level_heights(grid, Face(), Center()),   ("xf", "z"); attrib = Dict("long_name" => "height of u", "units" => "m"))
        defVar(dataset, "height_w", level_heights(grid, Center(), Face()),   ("x", "zf"); attrib = Dict("long_name" => "height of w", "units" => "m"))

        single(field) = Array{Float32}(Array(interior(field)))

        defVar(dataset, "u", single(u), ("xf", "y", "z"); attrib = Dict("long_name" => "wind component along x", "units" => "m s-1"))
        defVar(dataset, "v", single(v), ("x", "yf", "z"); attrib = Dict("long_name" => "wind component along y", "units" => "m s-1"))
        defVar(dataset, "w", single(w), ("x", "y", "zf"); attrib = Dict("long_name" => "vertical wind component", "units" => "m s-1"))
        defVar(dataset, "t", single(θ), ("x", "y", "z");  attrib = Dict("long_name" => "potential temperature", "units" => "K"))
        defVar(dataset, "p", single(p), ("x", "y", "z");  attrib = Dict("long_name" => "pressure", "units" => "Pa"))
    end

    return path
end
