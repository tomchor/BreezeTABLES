# Setup of the TABLES large-eddy simulation intercomparison cases with Breeze.
#
# The protocol is "Information on the Model Set-up and Data output" by Rauchöcker, Serafin and
# Stiperski (5 June 2025), stored in this repository as Rauchoecker_QJRMS_TABLES.pdf. Section
# numbers in the comments below refer to that document.
#
# Configuration: fully compressible dynamics with acoustic substepping, a terrain-following
# (Gal-Chen) vertical coordinate for the ridge cases, fifth-order WENO advection with no explicit
# subfilter-scale closure (implicit large-eddy simulation), and a neutral log-law wall model.
#
# TODO: once all cases work, consider adding an explicit subfilter-scale closure. The protocol
# requests subfilter-scale stresses, heat fluxes and dissipation (Section 3), which an implicit
# closure does not provide.

using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: znode
using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ, Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Breeze.TerrainFollowingDiscretization: TerrainFollowingVerticalDiscretization, LinearDecay,
                                             materialize_terrain!
using Printf
using Random

#####
##### Initial profiles (Section 2)
#####

"""
    MixedLayerSounding(FT; kwargs...)

Potential temperature of the flat and 100 m ridge cases: a mixed layer capped by an inversion,
with constant stratification above.

The protocol specifies a 4 K inversion on top of a 1000 m deep mixed layer but not its thickness.
The default of 50 m is read from Figure 1 of the protocol and should be replaced by the value in
the sounding files `input_sounding_teamx_*` once they are available.
"""
struct MixedLayerSounding{FT}
    mixed_layer_potential_temperature :: FT
    mixed_layer_depth :: FT
    inversion_strength :: FT
    inversion_thickness :: FT
    stratification :: FT
end

MixedLayerSounding(FT = Float64;
                   mixed_layer_potential_temperature = 300,
                   mixed_layer_depth = 1000,
                   inversion_strength = 4,
                   inversion_thickness = 50,
                   stratification = 0.003) =
    MixedLayerSounding{FT}(mixed_layer_potential_temperature, mixed_layer_depth,
                           inversion_strength, inversion_thickness, stratification)

@inline function (s::MixedLayerSounding)(z)
    z₁ = s.mixed_layer_depth
    z₂ = z₁ + s.inversion_thickness
    θ₁ = s.mixed_layer_potential_temperature
    θ₂ = θ₁ + s.inversion_strength
    θ_inversion = θ₁ + s.inversion_strength * (z - z₁) / s.inversion_thickness
    θ_above = θ₂ + s.stratification * (z - z₂)
    return ifelse(z < z₁, θ₁, ifelse(z < z₂, θ_inversion, θ_above))
end

"""
    StratifiedSounding(FT; kwargs...)

Potential temperature of the 1000 m ridge cases: constant stratification from `z = 0`.
"""
struct StratifiedSounding{FT}
    surface_potential_temperature :: FT
    stratification :: FT
end

StratifiedSounding(FT = Float64; surface_potential_temperature = 300, stratification = 0.003) =
    StratifiedSounding{FT}(surface_potential_temperature, stratification)

@inline (s::StratifiedSounding)(z) = s.surface_potential_temperature + s.stratification * z

"Wind in the x-direction that does not vary with height."
struct UniformWind{FT}
    speed :: FT
end

@inline (wind::UniformWind)(z) = wind.speed

"""
Wind in the x-direction that is zero below `base`, increases linearly to `speed` at `top`,
and is uniform above (forced convection over the 1000 m ridge).
"""
struct ShearLayerWind{FT}
    speed :: FT
    base :: FT
    top :: FT
end

@inline (wind::ShearLayerWind)(z) = wind.speed * clamp((z - wind.base) / (wind.top - wind.base), 0, 1)

#####
##### Cases (Section 2)
#####

const tables_case_names = (:flat_u00, :flat_u10, :ridge100_u00, :ridge100_u10, :ridge1000_u00, :ridge1000_u10)

"""
    tables_case(name, FT=Float64)

Return the parameters of the TABLES case `name`, one of `tables_case_names`: the ridge height,
the geostrophic wind, and the initial profiles of potential temperature and wind.
"""
function tables_case(name::Symbol, FT = Float64)
    mixed_layer = MixedLayerSounding(FT)
    stratified = StratifiedSounding(FT)
    calm = UniformWind{FT}(0)
    windy = UniformWind{FT}(10)
    sheared = ShearLayerWind{FT}(10, 1000, 1400)

    name === :flat_u00      && return (; name, ridge_height = 0,    geostrophic_wind = 0,  sounding = mixed_layer, wind = calm)
    name === :flat_u10      && return (; name, ridge_height = 0,    geostrophic_wind = 10, sounding = mixed_layer, wind = windy)
    name === :ridge100_u00  && return (; name, ridge_height = 100,  geostrophic_wind = 0,  sounding = mixed_layer, wind = calm)
    name === :ridge100_u10  && return (; name, ridge_height = 100,  geostrophic_wind = 10, sounding = mixed_layer, wind = windy)
    name === :ridge1000_u00 && return (; name, ridge_height = 1000, geostrophic_wind = 0,  sounding = stratified,  wind = calm)
    name === :ridge1000_u10 && return (; name, ridge_height = 1000, geostrophic_wind = 10, sounding = stratified,  wind = sheared)

    throw(ArgumentError("Unknown TABLES case $name. The cases are $tables_case_names."))
end

#####
##### Grid (Sections 1.1, 1.2 and 1.8)
#####

"""
    vertical_faces(; height=5000, spacing=10, uniform_height=3000, stretching=1.05, maximum_spacing=40)

Return the heights of the cell faces of the reference vertical coordinate: uniform `spacing` up to
`uniform_height`, then spacings that grow by the factor `stretching` per level up to
`maximum_spacing`. The stretched spacings are rescaled by a common factor so that the top face
is at `height`.

The defaults keep 10 m spacing up to the base of the damping layer at 3 km, which is also the top
of the common grid requested for the output (Section 5).
"""
function vertical_faces(; height = 5000, spacing = 10, uniform_height = 3000, stretching = 1.05, maximum_spacing = 40)
    uniform_levels = round(Int, uniform_height / spacing)
    uniform_levels * spacing ≈ uniform_height ||
        throw(ArgumentError("uniform_height = $uniform_height is not a multiple of spacing = $spacing"))

    faces = collect(range(0.0, uniform_height, length = uniform_levels + 1))
    remaining = height - uniform_height
    remaining > 0 || return faces

    spacings = Float64[]
    while sum(spacings) < remaining
        push!(spacings, min(spacing * stretching^(length(spacings) + 1), maximum_spacing))
    end
    spacings .*= remaining / sum(spacings)

    append!(faces, uniform_height .+ cumsum(spacings))
    faces[end] = height

    return faces
end

"""
    tables_grid(arch, FT, case; Nx=512, Ny=1024, Δx=20, z_faces=vertical_faces(), halo=(3, 3, 3))

Return the grid of a TABLES `case`. The domain is periodic in `x` and `y`. Cases with a ridge use a
terrain-following vertical coordinate with linear decay (Gal-Chen and Somerville 1975), for which
the vertical spacing is the reference spacing times `1 - h / height`. The ridge is
`h(x) = ridge_height / 2 * (1 - cos(2π x / Lx))`, with the valley at `x = 0`.

The default `halo` of three cells is what `WENO(order=5)` requires.
"""
function tables_grid(arch, FT, case; Nx = 512, Ny = 1024, Δx = 20, z_faces = vertical_faces(), halo = (3, 3, 3))
    Lx = Nx * Δx
    Ly = Ny * Δx
    Nz = length(z_faces) - 1
    horizontal = (; x = (0, Lx), y = (0, Ly), halo, topology = (Periodic, Periodic, Bounded))

    if case.ridge_height == 0
        return RectilinearGrid(arch, FT; size = (Nx, Ny, Nz), z = z_faces, horizontal...)
    end

    z = TerrainFollowingVerticalDiscretization(z_faces; formulation = LinearDecay())
    grid = RectilinearGrid(arch, FT; size = (Nx, Ny, Nz), z, horizontal...)

    h = case.ridge_height
    ridge(x, y) = h / 2 * (1 - cos(2π * x / Lx))
    materialize_terrain!(grid, ridge)

    return grid
end

#####
##### Wall model and surface heat flux (Section 1.4)
#####
##### Neutral log law applied locally at the first cell center, at height z₁ above the surface:
#####
#####   u★ = κ |U(z₁)| / log(z₁ / z₀),   τ = - ρ u★² U / |U|
#####
##### This is the definition of u★ used for the output variable `ust` (Section 3.1).
#####
##### TODO: replace the neutral log law with a stability-corrected Monin-Obukhov wall model
##### (Obukhov length from the prescribed heat flux) once the simulations are confirmed to work.
#####

@inline log_law_drag_coefficient(z₁, p) = (p.von_karman_constant / log(z₁ / p.roughness_length))^2

@inline function x_momentum_flux(i, j, grid, clock, fields, p)
    ρ = ℑxᶠᵃᵃ(i, j, 1, grid, fields.ρᵈ)
    u = @inbounds fields.u[i, j, 1]
    v = ℑxyᶠᶜᵃ(i, j, 1, grid, fields.v)
    z₁ = Δzᶠᶜᶜ(i, j, 1, grid) / 2
    return - ρ * log_law_drag_coefficient(z₁, p) * sqrt(u^2 + v^2) * u
end

@inline function y_momentum_flux(i, j, grid, clock, fields, p)
    ρ = ℑyᵃᶠᵃ(i, j, 1, grid, fields.ρᵈ)
    u = ℑxyᶜᶠᵃ(i, j, 1, grid, fields.u)
    v = @inbounds fields.v[i, j, 1]
    z₁ = Δzᶜᶠᶜ(i, j, 1, grid) / 2
    return - ρ * log_law_drag_coefficient(z₁, p) * sqrt(u^2 + v^2) * v
end

# The protocol prescribes the kinematic flux w′θ′. The flux of ρθ uses the density of the first cell.
@inline potential_temperature_flux(i, j, grid, clock, fields, p) =
    @inbounds fields.ρᵈ[i, j, 1] * p.kinematic_heat_flux

#####
##### Rayleigh damping toward the initial profiles (Section 1.3)
#####

@inline function rayleigh_damping_rate(z, p)
    s = clamp((z - p.damping_base) / (p.domain_height - p.damping_base), 0, 1)
    return p.damping_coefficient * (1 - cos(π * s)) / 2
end

@inline function u_damping(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, Face(), Center(), Center())
    return @inbounds - rayleigh_damping_rate(z, p) * (fields.u[i, j, k] - p.wind(z))
end

@inline function v_damping(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, Center(), Face(), Center())
    return @inbounds - rayleigh_damping_rate(z, p) * fields.v[i, j, k]
end

@inline function w_damping(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, Center(), Center(), Face())
    return @inbounds - rayleigh_damping_rate(z, p) * fields.w[i, j, k]
end

@inline function θ_damping(i, j, k, grid, clock, fields, p)
    z = znode(i, j, k, grid, Center(), Center(), Center())
    return @inbounds - rayleigh_damping_rate(z, p) * (fields.θ[i, j, k] - p.sounding(z))
end

#####
##### Model
#####

"""
    tables_model(grid, case; kwargs...)

Return an `AtmosphereModel` for the TABLES `case` on `grid`. The model has no explicit
subfilter-scale closure: the dissipation of the advection scheme takes its place.

Keyword arguments
=================

- `advection`: advection scheme for momentum and scalars. Default: `WENO(order=5)`.
- `time_discretization`: default `SplitExplicitTimeDiscretization()`, Wicker-Skamarock RK3 with
  acoustic substeps chosen from an acoustic Courant number of 0.5.
- `roughness_length` (m), `kinematic_heat_flux` (K m s⁻¹), `base_pressure` (Pa): Section 1.4.
- `coriolis_parameter` (s⁻¹): Section 1.7.
- `damping_base` (m), `damping_coefficient` (s⁻¹): Section 1.3.
- `timestepper_kwargs`: forwarded to `AtmosphereModel`, for example
  `(; rounding = StochasticRounding())` for single precision.
"""
function tables_model(grid, case;
                      advection = WENO(order = 5),
                      time_discretization = SplitExplicitTimeDiscretization(),
                      roughness_length = 0.1,
                      kinematic_heat_flux = 0.12,
                      base_pressure = 1e5,
                      coriolis_parameter = 1e-4,
                      damping_base = 3000,
                      damping_coefficient = 1/300,
                      timestepper_kwargs = NamedTuple())

    FT = eltype(grid)
    sounding = case.sounding
    wind = case.wind

    # The reference state is the initial sounding in discrete hydrostatic balance.
    reference_potential_temperature(z) = sounding(z)
    dynamics = CompressibleDynamics(time_discretization;
                                    base_pressure = FT(base_pressure),
                                    standard_pressure = FT(1e5),
                                    reference_potential_temperature)

    # Coriolis force and geostrophic pressure gradient (Section 1.7)
    coriolis = FPlane(FT; f = coriolis_parameter)
    geostrophic_wind = FT(case.geostrophic_wind)
    uᵍ(z) = geostrophic_wind
    vᵍ(z) = zero(geostrophic_wind)
    geostrophic = geostrophic_forcings(uᵍ, vᵍ)

    # Rayleigh damping (Section 1.3)
    damping = (; damping_coefficient = FT(damping_coefficient),
                 damping_base = FT(damping_base),
                 domain_height = FT(grid.Lz),
                 sounding, wind)

    damp(func) = Forcing(func; discrete_form = true, parameters = damping)

    forcing = (u = (geostrophic.u, damp(u_damping)),
               v = (geostrophic.v, damp(v_damping)),
               w = damp(w_damping),
               θ = damp(θ_damping))

    # Surface fluxes (Section 1.4)
    wall = (; von_karman_constant = FT(0.4), roughness_length = FT(roughness_length))
    heating = (; kinematic_heat_flux = FT(kinematic_heat_flux))

    ρu_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(x_momentum_flux; discrete_form = true, parameters = wall))
    ρv_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(y_momentum_flux; discrete_form = true, parameters = wall))
    ρθ_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(potential_temperature_flux; discrete_form = true, parameters = heating))

    boundary_conditions = (ρu = ρu_bcs, ρv = ρv_bcs, ρθ = ρθ_bcs)

    return AtmosphereModel(grid; dynamics, coriolis, advection, forcing,
                           boundary_conditions, timestepper_kwargs)
end

#####
##### Initial condition (Sections 1.5 and 2)
#####

"""
    set_tables_initial_condition!(model, case; noise_amplitude=0.25, seed=2026)

Set the initial sounding of `case` plus potential temperature noise drawn uniformly from
`[-noise_amplitude, noise_amplitude]` at every cell (Section 1.5).

The density is set to `ρᵣ θᵣ / θ`, where `ρᵣ` and `θᵣ` are the reference density and the sounding.
This leaves `ρθ`, and therefore the pressure, equal to the hydrostatically balanced reference,
so the noise produces buoyancy but no acoustic transient.
"""
function set_tables_initial_condition!(model, case; noise_amplitude = 0.25, seed = 2026)
    grid = model.grid
    arch = architecture(grid)
    FT = eltype(grid)
    sounding = case.sounding
    wind = case.wind

    # ρ temporarily holds the sounding θᵣ
    ρ = CenterField(grid)
    set!(ρ, (x, y, z) -> sounding(z))

    noise = rand(Xoshiro(seed), FT, size(grid)...)
    noise = on_architecture(arch, noise)
    θ = CenterField(grid)
    interior(θ) .= interior(ρ) .+ FT(noise_amplitude) .* (2 .* noise .- 1)

    ρᵣ = model.dynamics.reference_state.density
    interior(ρ) .= interior(ρᵣ) .* interior(ρ) ./ interior(θ)

    set!(model; ρ, θ, u = (x, y, z) -> wind(z), v = 0, w = 0, enforce_mass_conservation = false)

    return nothing
end

"""
    maximum_pressure_imbalance(model)

Return `max|p - pᵣ|` in Pa. This is at round-off level after `set_tables_initial_condition!`.
"""
function maximum_pressure_imbalance(model)
    p = model.dynamics.pressure
    pᵣ = model.dynamics.reference_state.pressure
    return maximum(abs, interior(p) .- interior(pᵣ))
end

#####
##### Simulation
#####

"""
    tables_simulation(model; Δt=0.2, stop_time=3hours, progress_interval=500)

Return a `Simulation` with a fixed time step (Section 1.6) and a progress message.
"""
function tables_simulation(model; Δt = 0.2, stop_time = 3hours, progress_interval = 500)
    simulation = Simulation(model; Δt, stop_time, minimum_relative_step = 1e-6)
    Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)

    wall_clock = Ref(time_ns())

    function progress(sim)
        u, v, w = sim.model.velocities
        elapsed = 1e-9 * (time_ns() - wall_clock[])
        @info @sprintf("iteration %d, time %s, wall time %s, max|u| %.2f, max|v| %.2f, max|w| %.2f m/s",
                       iteration(sim), prettytime(sim), prettytime(elapsed),
                       maximum(abs, u), maximum(abs, v), maximum(abs, w))
        return nothing
    end

    add_callback!(simulation, progress, IterationInterval(progress_interval))

    return simulation
end
