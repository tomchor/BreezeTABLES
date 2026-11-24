# Animate the x-z output of a TABLES case written by run_case.jl.
#
#     julia --project animate.jl output_directory case_name
#
# Reads `case_name_slices.jld2` (instantaneous slices) and `case_name_sections.jld2` (averages along y)
# from `output_directory` and writes `case_name.mp4` and `case_name_final.png` there. The two files
# must have been saved at the same times, which is the default of run_case.jl.
#
# Top panel: vertical velocity in the slice. Bottom panel: wind across the ridge, averaged along
# the ridge. Lines are isentropes. Fields are interpolated from model levels to constant heights.

using Breeze # defines the terrain-following grid stored in the output files
using Oceananigans
using Oceananigans: instantiated_location
using Oceananigans.Grids: znode
using CairoMakie
using Printf
using Statistics

output_directory = ARGS[1]
case_name = ARGS[2]

slices_file = joinpath(output_directory, case_name * "_slices.jld2")
sections_file = joinpath(output_directory, case_name * "_sections.jld2")

w = FieldTimeSeries(slices_file, "w")
θ = FieldTimeSeries(slices_file, "θ")
u_mean = FieldTimeSeries(sections_file, "u")
θ_mean = FieldTimeSeries(sections_file, "θ")

times = w.times
times == u_mean.times || error("The slices and the sections were not saved at the same times.")

grid = w.grid
Lx, Lz = grid.Lx, grid.Lz

#####
##### Interpolation from model levels to constant heights
#####

# The terrain does not vary along y, so the heights are evaluated at j = 1.
node_heights(grid, ℓx, ℓz, Nx, Nz) = [znode(i, 1, k, grid, ℓx, Center(), ℓz) for i in 1:Nx, k in 1:Nz]
surface_heights(grid, ℓx, Nx) = [znode(i, 1, 1, grid, ℓx, Center(), Face()) for i in 1:Nx]

"""
    at_constant_heights(data, heights, surface, levels; fill_depth = 0)

Interpolate each column of `data`, given at `heights[i, k]`, linearly to `levels`. Levels below the
lowest node take the value of that node, down to `fill_depth` below `surface[i]`. Levels further
below are `NaN`.
"""
function at_constant_heights(data, heights, surface, levels; fill_depth = 0)
    Nx, Nz = size(data)
    interpolated = fill(NaN32, Nx, length(levels))

    for i in 1:Nx, (n, z) in enumerate(levels)
        z < surface[i] - fill_depth && continue

        k = searchsortedlast(view(heights, i, :), z)
        if k == 0
            interpolated[i, n] = data[i, 1]
        elseif k == Nz
            interpolated[i, n] = data[i, Nz]
        else
            weight = (z - heights[i, k]) / (heights[i, k+1] - heights[i, k])
            interpolated[i, n] = (1 - weight) * data[i, k] + weight * data[i, k+1]
        end
    end

    return interpolated
end

"Return the x nodes of `series`, the terrain height there, and every saved field at constant heights."
function regrid(series, levels)
    ℓx, ℓy, ℓz = instantiated_location(series[1])
    Nx, _, Nz = size(interior(series[1]))
    heights = node_heights(grid, ℓx, ℓz, Nx, Nz)
    surface = surface_heights(grid, ℓx, Nx)

    # Each column is drawn as a rectangle, which leaves gaps next to sloping terrain. Filling a
    # little below the surface closes them; the terrain drawn on top hides the filled part.
    fill_depth = maximum(abs, diff(surface)) / 2 + step(levels)

    frames = [at_constant_heights(dropdims(interior(series[n]), dims = 2), heights, surface, levels; fill_depth)
              for n in 1:length(series.times)]
    return xnodes(grid, ℓx), surface, frames
end

plot_top = min(Lz, 3000) # the damping layer starts at 3 km
levels = range(10, plot_top, step = 20)

x_w, _, w_frames = regrid(w, levels)
x_θ, surface, θ_frames = regrid(θ, levels)
x_u, _, u_frames = regrid(u_mean, levels)
_, _, θ_mean_frames = regrid(θ_mean, levels)

#####
##### Figure
#####

km = 1e-3
finite(frames) = filter(isfinite, reduce(vcat, vec.(frames)))
w_limit = max(0.1, round(quantile(abs.(finite(w_frames)), 0.999), sigdigits = 1))
u_limit = max(0.1, round(quantile(abs.(finite(u_frames)), 0.999), sigdigits = 1))
isentropes = 280:0.5:330

n = Observable(1)

Δx = round(Int, grid.Lx / grid.Nx)
Δz = round(Int, first(diff(znodes(grid, Face())))) # reference spacing of the lowest level

title = @lift @sprintf("TABLES %s with %d m × %d m cells, t = %s", case_name, Δx, Δz, prettytime(times[$n]))

fig = Figure(size = (1100, 760), fontsize = 15)
Label(fig[0, 1:2], title, fontsize = 17, tellwidth = false)

panels = ((1, "Vertical velocity in an x-z slice", x_w, w_frames, θ_frames, w_limit, "w (m s⁻¹)"),
          (2, "Wind across the ridge, averaged along the ridge", x_u, u_frames, θ_mean_frames, u_limit, "u (m s⁻¹)"))

for (row, panel_title, x, frames, θ_contours, limit, label) in panels
    ax = Axis(fig[row, 1]; title = panel_title, ylabel = "z (km)", aspect = DataAspect(),
              limits = ((0, Lx * km), (0, plot_top * km)))
    row == 2 ? (ax.xlabel = "x (km)") : hidexdecorations!(ax, grid = false)

    hm = heatmap!(ax, x .* km, levels .* km, @lift(frames[$n]);
                  colormap = :balance, colorrange = (-limit, limit), nan_color = :transparent)
    contour!(ax, x_θ .* km, levels .* km, @lift(θ_contours[$n]);
             levels = isentropes, color = (:black, 0.55), linewidth = 0.6)
    band!(ax, x_θ .* km, zero(surface), surface .* km; color = :gray35)

    Colorbar(fig[row, 2], hm; label)
end

#####
##### Animation
#####

movie_file = joinpath(output_directory, case_name * ".mp4")
record(fig, movie_file, eachindex(times); framerate = 15) do frame
    n[] = frame
end

n[] = length(times)
image_file = joinpath(output_directory, case_name * "_final.png")
save(image_file, fig)

@info "Wrote $movie_file and $image_file"
