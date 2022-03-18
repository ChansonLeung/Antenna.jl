import PlotlyJS
using LinearAlgebra
using Antenna
using Interpolations

function ploting_point(p::Vector{anten_point})
    expand_point = (vec_point) -> begin
        ([i.p.x for i = vec_point],
            [i.p.y for i = vec_point],
            [i.p.z for i = vec_point])
    end
    x, y, z = expand_point(p)

    max = maximum(sqrt.(x .^ 2 + y .^ 2 + z .^ 2))

    PlotlyJS.plot(PlotlyJS.scatter(
            x = x, y = y, z = z, text = 1:size(x, 1),
            type = "scatter3d",
            marker = PlotlyJS.attr(
                size = 2
            ),
            mode = "markers",
            hovertemplate = "x:%{x:.3f} <br>y:%{y:.3f} <br>x:%{z:.3f} <br>i:%{text} <extra></extra> ",
        ),
        PlotlyJS.Layout(
            scene_xaxis_range = [-max, max],
            scene_yaxis_range = [-max, max],
            scene_zaxis_range = [-max, max])
    )
end
function ploting_point(p::Vector{Vector{Float64}})
    expand_point = (vec_point) -> begin
        ([i[1] for i = vec_point],
            [i[2] for i = vec_point],
            [i[3] for i = vec_point])
    end
    x, y, z = expand_point(p)
    max = maximum(sqrt.(x .^ 2 + y .^ 2 + z .^ 2))
    PlotlyJS.plot(PlotlyJS.scatter(x = x, y = y, z = z,
        type = "scatter3d",
        marker = PlotlyJS.attr(
            size = 2
        ),
        PlotlyJS.Layout(
            xaxis_range = [-max, max],
            yaxis_range = [-max, max],
            zaxis_range = [-max, max]),
        mode = "markers"))
end
# ploting_point(point)

function ploting_pattern(pattern::anten_pattern; min = -40)
    r = directivity(pattern)

    r_log_raw = 10log10.(r)
    # get maximum U and the direction
    max_U = maximum(r_log_raw)
    max_theta = θ_grid[r_log_raw.==max_U] .|>
                rad2deg .|>
                x -> round(x, digits = 2)
    max_phi = ϕ_grid[r_log_raw.==max_U] .|>
              rad2deg .|>
              x -> round(x, digits = 2)

    r_log_limmin = 10log10.(r)
    r_log_limmin[r_log_limmin.<min] .= min
    r_plot = 10log10.(r)
    r_plot[r_plot.<min] .= min
    r_plot .-= min

    θ = [θ for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    ϕ = [ϕ for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    tulple_matrix = sph2cart.(θ, ϕ, r_plot)
    x = [i.x for i = tulple_matrix]
    y = [i.y for i = tulple_matrix]
    z = [i.z for i = tulple_matrix]

    max = maximum(r_plot)
    PlotlyJS.plot(
        PlotlyJS.surface(
            x = x, y = y, z = z,
            # customdata = [[θ,ϕ] for (θ,ϕ) in zip(rad2deg.(θ), rad2deg.(ϕ))],
            text = map(
                p -> "θ:$(round(p[1], digits=2))    ϕ:$(round(p[2], digits=2))    r:$(round(p[3], digits=2))",
                zip(
                    rad2deg.(θ),
                    rad2deg.(ϕ),
                    r_log_raw
                )
            ),
            surfacecolor = r_log_limmin,
            colorscale = "Jet",
            # hovertemplate = "x:%{x},y:%{y},z:%{z}, <br>θ:%{customdata} <br>ϕ:%{customdata[1]} ",
            # hovertemplate = "x:%{x},y:%{y},z:%{z}, <br>θ:%{customdata}",
        ),
        PlotlyJS.Layout(
            title = PlotlyJS.attr(
                text = "最大方向性系数: $(round(max_U, digits=2)) dB in θ=$(max_theta)°,ϕ=$(max_phi)°"
            ),
            scene = PlotlyJS.attr(
                xaxis = PlotlyJS.attr(
                    showgrid = false,
                    showticklabels = false,
                    showbackground = false,
                    showaxeslabels = false,
                    # title= PlotlyJS.attr(
                    #     text=""
                    # )
                ),
                yaxis = PlotlyJS.attr(
                    showgrid = false,
                    showticklabels = false,
                    showbackground = false,
                    # title= PlotlyJS.attr(
                    #     text=""
                    # )
                ),
                zaxis = PlotlyJS.attr(
                    showgrid = false,
                    # showbackground = false,
                    showticklabels = false,
                    title = PlotlyJS.attr(
                        text = ""
                    )
                ),
            ),
            scene_xaxis_range = [-max, max],
            scene_yaxis_range = [-max, max],
            scene_zaxis_range = [-max, max]
        )
    )

end
function ploting_pattern_2D(pattern::anten_pattern; θ = deg2rad.(-179:179), ϕ = 0)
    # θ ∈ [0 180]
    # ϕ ∈ [-180 180]
    sign = x-> x>=0 ? 1 : -1
    item1 = -sign.(mod2pi.(θ) .- π)
    θ_ = mod.(θ.*item1 , π)
    ϕ_ = item1 .* ϕ .|> mod2pi
    r = 1 / 2(120pi) * (abs.(pattern.θ.(θ_, ϕ_)).^2 .+ abs.(pattern.ϕ.(θ_, ϕ_)).^2)
    r = 4pi * r / radiated_power(pattern)
    r = 10log10.(r)
    PlotlyJS.plot(rad2deg.(θ), r)
end
# ploting_pattern(global_pattern)

function plot_E_vec_field(pattern::anten_pattern; mode = "1vec")
    # get point from default θ and ϕ
    θ = [θ for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    ϕ = [ϕ for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    point_tuple_matrix = sph2cart.(θ, ϕ, 1)
    x = [i.x for i = vec(point_tuple_matrix)]
    y = [i.y for i = vec(point_tuple_matrix)]
    z = [i.z for i = vec(point_tuple_matrix)]

    vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
    vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0]
    # get pattern  E vector

    #convert E vector to 1D vector

    if mode == "2vec"
        vec_vec_θ = [vec_θ(θ, ϕ) * pattern.θ(θ, ϕ) for (θ, ϕ) = vec(zip(θ_grid, ϕ_grid))]
        vec_vec_ϕ = [vec_ϕ(θ, ϕ) * pattern.ϕ(θ, ϕ) for (θ, ϕ) = vec(zip(θ_grid, ϕ_grid))]
        u = [i[1] for i = vec_vec_θ]
        v = [i[2] for i = vec_vec_θ]
        w = [i[3] for i = vec_vec_θ]
        u2 = [i[1] for i = vec_vec_ϕ]
        v2 = [i[2] for i = vec_vec_ϕ]
        w2 = [i[3] for i = vec_vec_ϕ]
        x = [x; x]
        y = [y; y]
        z = [z; z]
        u = [u; u2]
        v = [v; v2]
        w = [w; w2]
    else
        vec_vec_globe = [
            vec_θ(θ, ϕ) * abs(pattern.θ(θ, ϕ)) .+ vec_ϕ(θ, ϕ) * abs(pattern.ϕ(θ, ϕ))
            for (θ, ϕ) = zip(θ_grid, ϕ_grid)
        ]

        # vec_vec_globe_1 = [
        #     vec_θ(θ, ϕ) * pattern.θ(θ, ϕ) .+ vec_ϕ(θ, ϕ) * pattern.ϕ(θ, ϕ)
        #     for (θ,ϕ) = zip(θ_grid, ϕ_grid)
        # ]

        # @show vec_vec_globe_1 = [
        #     [vec_θ(θ, ϕ) , vec_ϕ(θ, ϕ)]
        #     for (θ,ϕ) = zip(deg2rad(45),deg2rad(-90))
        # ]
        # @show minimum(norm.(vec_vec_globe))
        # @show minimum(norm.(vec_vec_globe_1))
        # @show minimum((vec_vec_globe_1))
        # vec_vec_globe = [[1,1,1] for (θ,ϕ) = zip(vec(θ_grid), vec(ϕ_grid))]
        u = [i[1] for i = vec(vec_vec_globe)]
        v = [i[2] for i = vec(vec_vec_globe)]
        w = [i[3] for i = vec(vec_vec_globe)]
    end
    PlotlyJS.plot(
        cone(
            x = x,
            y = y,
            z = z,
            u = u,
            v = v,
            w = w,
            sizemode = "scaled",
            sizeref = 1
        )
    )

    # PlotlyJS.plot(PlotlyJS.scatter(x=x,y=y,z=z, 
    #             type="scatter3d",
    #              marker=PlotlyJS.attr(
    #                  size=2
    #             ),
    #              mode="markers"))
end

export
    ploting_point,
    ploting_pattern,
    plot_E_vec_field,
    plot_E_vec_field_test,
    ploting_pattern_2D