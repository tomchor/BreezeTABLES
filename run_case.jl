# Run one TABLES case.
#
#     julia --project run_case.jl flat_u00
#
# The case is one of flat_u00, flat_u10, ridge100_u00, ridge100_u10, ridge1000_u00, ridge1000_u10.
# The defaults follow the protocol at full size (512 x 1024 columns at 20 m, time step 0.2 s, 3 h).
# Environment variables override them, mainly to run reduced configurations for testing:
#
#     TABLES_ARCH            cpu or gpu (default: gpu)
#     TABLES_FLOAT_TYPE      Float64 or Float32 (default: Float64)
#     TABLES_NX, TABLES_NY   number of columns in x and y (default: 512, 1024)
#     TABLES_DX              horizontal spacing in m (default: 20)
#     TABLES_DZ              vertical spacing in m below TABLES_UNIFORM_HEIGHT (default: 10)
#     TABLES_UNIFORM_HEIGHT  height in m below which the vertical spacing is uniform (default: 3000)
#     TABLES_MAX_DZ          largest vertical spacing in m (default: 40)
#     TABLES_DT              time step in s (default: 0.2)
#     TABLES_STOP_TIME       duration in s (default: 10800)
#     TABLES_OUTPUT_DIR      output directory (default: output)
#     TABLES_PICKUP          true to restart from the latest checkpoint (default: false)

include(joinpath(@__DIR__, "tables.jl"))

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

simulation = tables_simulation(model; Δt, stop_time)

#####
##### Output
#####
##### This is monitoring output for checking that a case runs correctly. The statistics, snapshots
##### and spectra requested in Sections 3 to 5 of the protocol are not implemented yet.
#####

pickup = get(ENV, "TABLES_PICKUP", "false") == "true"

iterations(interval) = max(1, round(Int, interval / Δt))

u, v, w = model.velocities
θ = model.formulation.potential_temperature

outputs = (; u, v, w, θ, ww = w^2, wθ = w * θ)

# Averages along the ridge, on model levels, every 10 minutes
sections = NamedTuple(name => Average(@at((Center, Center, Center), outputs[name]), dims = 2)
                      for name in keys(outputs))

simulation.output_writers[:sections] = JLD2Writer(model, sections;
                                                  dir = output_dir,
                                                  filename = "$(case_name)_sections",
                                                  schedule = IterationInterval(iterations(10minutes)),
                                                  overwrite_files = !pickup)

# Checkpoints every 30 minutes. The one at 2.5 h is the start of the statistics window.
simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                        dir = output_dir,
                                                        prefix = "$(case_name)_checkpoint",
                                                        schedule = IterationInterval(iterations(30minutes)))

run!(simulation; pickup)
