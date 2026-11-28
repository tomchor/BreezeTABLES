# Animate the slices and sections of a TABLES case written by run_case.jl.
#
#     julia --project animate.jl output_directory case_name
#
# Reads `case_name_slices.jld2` (instantaneous x-z slices), `case_name_xy_slices.jld2` (instantaneous
# x-y slices) and `case_name_sections.jld2` (averages along y) from `output_directory` and writes
# `case_name.mp4` and `case_name_final.png` there. The three files must have been saved at the same
# times, which is the default of run_case.jl.
#
# Left panel: vertical velocity in the x-y slice, which is on a model level. Top right panel: vertical
# velocity in the x-z slice. Bottom right panel: wind across the ridge, averaged along the ridge.
# Solid lines are isentropes. Dashed lines show where each slice cuts the other. The x-z fields are
# interpolated from model levels to constant heights.

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
xy_slices_file = joinpath(output_directory, case_name * "_xy_slices.jld2")
sections_file = joinpath(output_directory, case_name * "_sections.jld2")

w = FieldTimeSeries(slices_file, "w")
θ = FieldTimeSeries(slices_file, "θ")
w_xy = FieldTimeSeries(xy_slices_file, "w")
u_mean = FieldTimeSeries(sections_file, "u")
θ_mean = FieldTimeSeries(sections_file, "θ")

times = w.times
times == w_xy.times || error("The x-z and the x-y slices were not saved at the same times.")
times == u_mean.times || error("The slices and the sections were not saved at the same times.")

grid = w.grid
Lx, Ly, Lz = grid.Lx, grid.Ly, grid.Lz

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
##### The x-y slice, which stays on its model level
#####

x_xy = xnodes(grid, Center())
y_xy = ynodes(grid, Center())
xy_frames = [interior(w_xy[n], :, :, 1) for n in eachindex(times)]

# Where each slice lies in the other
y_slice = y_xy[only(w.indices[2])]
level_heights = [znode(i, 1, only(w_xy.indices[3]), grid, Center(), Center(), Face()) for i in 1:grid.Nx]

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

xy_title = if maximum(level_heights) - minimum(level_heights) < 1
    @sprintf("Vertical velocity at z = %.0f m", first(level_heights))
else
    @sprintf("Vertical velocity on the model level\nat z = %.0f to %.0f m", extrema(level_heights)...)
end

"Size of the largest panel that fits in `width × height` with the aspect ratio of a `Lx × Ly` domain."
panel_size(Lx, Ly, width, height) = round.(Int, (min(width, height * Lx / Ly), min(height, width * Ly / Lx)))

# Every panel keeps the aspect ratio of the data. With the 10.24 km × 20.48 km of the protocol, the
# sizes give a figure of about 1280 × 720.
xy_width, xy_height = panel_size(Lx, Ly, 400, 570)
xz_width, xz_height = panel_size(Lx, plot_top, 1040 - xy_width, 240)

fig = Figure(fontsize = 15)
Label(fig[0, 1:2], title, fontsize = 17, tellwidth = false)

left = fig[1, 1] = GridLayout()
right = fig[1, 2] = GridLayout()
dashes = (color = :black, linestyle = :dash, linewidth = 1)

ax_xy = Axis(left[1, 1]; title = xy_title, xlabel = "x (km)", ylabel = "y (km)",
             width = xy_width, height = xy_height, limits = ((0, Lx * km), (0, Ly * km)))

# The colours are those of the x-z slice
heatmap!(ax_xy, x_xy .* km, y_xy .* km, @lift(xy_frames[$n]); colormap = :balance, colorrange = (-w_limit, w_limit))
hlines!(ax_xy, y_slice * km; dashes...)

panels = ((1, "Vertical velocity in an x-z slice", x_w, w_frames, θ_frames, w_limit, "w (m s⁻¹)"),
          (2, "Wind across the ridge, averaged along the ridge", x_u, u_frames, θ_mean_frames, u_limit, "u (m s⁻¹)"))

for (row, panel_title, x, frames, θ_contours, limit, label) in panels
    ax = Axis(right[row, 1]; title = panel_title, ylabel = "z (km)",
              width = xz_width, height = xz_height, limits = ((0, Lx * km), (0, plot_top * km)))
    row == 2 ? (ax.xlabel = "x (km)") : hidexdecorations!(ax, grid = false)

    hm = heatmap!(ax, x .* km, levels .* km, @lift(frames[$n]);
                  colormap = :balance, colorrange = (-limit, limit), nan_color = :transparent)
    contour!(ax, x_θ .* km, levels .* km, @lift(θ_contours[$n]);
             levels = isentropes, color = (:black, 0.55), linewidth = 0.6)
    band!(ax, x_θ .* km, zero(surface), surface .* km; color = :gray35)
    row == 1 && lines!(ax, x_xy .* km, level_heights .* km; dashes...)

    Colorbar(right[row, 2], hm; label)
end

resize_to_layout!(fig)

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
