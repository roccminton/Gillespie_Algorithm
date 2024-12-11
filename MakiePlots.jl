#Functions for Plotting using Makie

using CairoMakie

f = Figure()
ax = Axis(f[1, 1], xgridvisible = false, ygridvisible = false)
ylims!(ax, -1, 2)
bracket!(ax, 1, 0, 3, 0, text = "Curly", style = :curly)
bracket!(ax, 2, 1, 4, 1, text = "Square", style = :square)

f
