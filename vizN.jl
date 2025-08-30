"""
## Example

```
N0 = N1k(0)
N1 = N1k(1)
@variables u[1:2]
N0_f = [eval(build_function(N0[i], u)[1]) for i in eachindex(N0)]
f0 = plot_nedelec_functions_on_triangle(N0_f, density = 15)
N1_f = [eval(build_function(N1[i], u)[1]) for i in eachindex(N1)]
f1 = plot_nedelec_functions_on_triangle(N1_f, density = 12)
``` 
"""
function plot_nedelec_functions_on_triangle(functions; max_columns = 3, lengthscale = 0.05, density = 25)

    fig = Figure()

    points_triang = [Point2f(0.0), Point2f(0.0, 1.0), Point2f(1.0, 0.0)]

    n_points = density
    s = LinRange(0, 1, n_points)
    t = LinRange(0, 1, n_points)

    points = vec([(si, ti) for si in s, ti in t if (si + ti <= 1.0)])

    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]

    colors = [RGBf(rand(3)...) for _ in eachindex(functions)]
    cplot = 0
    i = 1
    for (idx, func) in enumerate(functions)
        idplot = mod(idx-1, max_columns)+1
        if idplot == 1
            cplot = cplot + 1
        end
        ax = Axis(fig[cplot, idplot], subtitle = Latexify.LaTeXString("\\phi_{$i}"), subtitlesize=18)
        i += 1
        poly!(ax, points_triang, color=:white, strokecolor=:black, strokewidth=2)

        ni = map(x -> func(x), points)

        u_coords = [p[1] for p in ni]
        v_coords = [p[2] for p in ni]

        arrows2d!(ax, x_coords, y_coords, u_coords, v_coords,
                color=colors[idx], lengthscale = lengthscale, label="N$idx")
    end

    return fig
end

