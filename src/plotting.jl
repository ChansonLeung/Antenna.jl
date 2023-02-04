import PlotlyJS
using LinearAlgebra
using Antenna
using Interpolations
using Plots
using Lazy:@>


@recipe f(::Type{anten_point}, point::anten_point) = point.pattern

@recipe f(::Type{Vector{anten_point}}, points::Vector{anten_point}) = getfield.(points, :p)

@recipe function f(pattern::anten_pattern, component=:all)
    xlabel --> "u"
    ylabel --> "v"

    uv_grid =[(u,v) for u in LinRange(-1,1,361), v in LinRange(-1,1,361)] 
    u = @> uv_grid getindex.(1)
    v = @> uv_grid getindex.(2)
    w(u,v) =  u^2+v^2 <=1 ? sqrt(1-u^2-v^2) : 0.
    w = w.(u,v)


    θ = @. acos(w/sqrt(u^2+v^2+w^2))
    ϕ = @. atan(v,u)
    r = @> directivity(pattern, component).(θ, ϕ)  .|> x->10log10(x)

    @. r[u^2+v^2>1] = NaN

    clims --> (-40, maximum(r))
    @series begin
        seriestype := :heatmap 
        (u[:,1], v[1,:], r)
    end
    # @series begin
    #     seriestype := :contour 
    #     (u[:,1],v[1,:],r)
    # end

end

@recipe function f(points::Vector{Vector{Float64}};)
    seriestype-->:scatter3d
    markersize-->1
    (getindex.(points,1),
    getindex.(points,2),
    getindex.(points,3))
end

function plot_point(p::Vector{Vector{Float64}}; ret_trace=false)
    expand_point = (vec_point) -> begin
        ([i[1] for i = vec_point],
            [i[2] for i = vec_point],
            [i[3] for i = vec_point])
    end
    x, y, z = expand_point(p)
    max = maximum(sqrt.(x .^ 2 + y .^ 2 + z .^ 2))
    trace = PlotlyJS.scatter(
            x = x, y = y, z = z, text = 1:size(x, 1),
            type = "scatter3d",
            marker = PlotlyJS.attr(
                size = 2
            ),
            mode = "markers",
            hovertemplate = "x:%{x:.3f} <br>y:%{y:.3f} <br>x:%{z:.3f} <br>i:%{text} <extra></extra> ",
            showlegend = false
        )
    layout_args = Dict(
        :scene_xaxis_range => [-max, max],
        :scene_yaxis_range => [-max, max],
        :scene_zaxis_range => [-max, max],
    )
    if ret_trace
        trace, layout_args
    else
        fig = PlotlyJS.plot(trace)
        PlotlyJS.relayout!(fig; layout_args...)
        display(fig)
    end
end
function plot_point(point::Vector{anten_point}; ret_trace = false)
    # expand tuple to vector
    vec_vec_point = map(point) do i  
        [i.p...] 
    end
    plot_point(vec_vec_point, ret_trace = ret_trace)
end

function plot_pattern(pattern::anten_pattern; min = -20, θ = θ_default, ϕ = ϕ_default, ret_trace=false)
    θ_grid = [θ for θ in θ, ϕ in ϕ]
    ϕ_grid = [ϕ for θ in θ, ϕ in ϕ]

    r = directivity(pattern).(θ_grid, ϕ_grid)
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
    x = [i[1] for i = tulple_matrix]
    y = [i[2] for i = tulple_matrix]
    z = [i[3] for i = tulple_matrix]

    max = maximum(r_plot)
    trace = PlotlyJS.surface(
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
        )
    layout_args = Dict(
                :title => PlotlyJS.attr(
                    text = "最大方向性系数: $(round(max_U, digits=2)) dB in θ=$(max_theta)°,ϕ=$(max_phi)°"
                ),
                :scene => PlotlyJS.attr(
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
                :scene_xaxis_range => [-max, max],
                :scene_yaxis_range => [-max, max],
                :scene_zaxis_range => [-max, max]
            )
    if ret_trace
        trace,layout_args
    else
        PlotlyJS.plot(
            trace, PlotlyJS.Layout(layout_args)
        )
    end

end
function plot_pattern_2D(pattern::anten_pattern; θ = rad2deg.([-reverse(θ_default);θ_default]), ϕ = 0, ret_trace=false)
    θ_ = [mod_angle_deg(θ,ϕ)[1] |>deg2rad for θ in θ]
    ϕ_ = [mod_angle_deg(θ,ϕ)[2] |>deg2rad for θ in θ]
    r = directivity(pattern).(θ_,ϕ_) |> x->10log10.(x)
    if ret_trace 
        PlotlyJS.scatter(x=θ,y=r, mode="lines", showlegend = false)
    else
        PlotlyJS.plot(θ, r)
    end

end
function plot_pattern_all(pattern::anten_pattern; θ = rad2deg.([-reverse(θ_default);θ_default]), ϕ = 0, point = false)
    trace_line = plot_pattern_2D(pattern, θ=θ, ϕ=ϕ, ret_trace=true)
    trace_3d, trace_3d_layout = plot_pattern(pattern, ret_trace=true)
    fig = empty
    if point isa Bool
        fig = PlotlyJS.make_subplots(
            rows=2, cols=2,
            specs=[
                PlotlyJS.Spec(kind="xy", colspan=2) missing
                PlotlyJS.Spec(kind="scene", colspan=1) missing 
            ]
        )
    else
        fig = PlotlyJS.make_subplots(
            rows=2, cols=2,
            specs=[
                PlotlyJS.Spec(kind="xy", colspan=2) missing
                PlotlyJS.Spec(kind="scene") PlotlyJS.Spec(kind="scene")
            ]
        )
        trace_point, trace_point_layout = plot_point(point, ret_trace=true)
        PlotlyJS.add_trace!(fig, trace_point, row=2,col=2)
        PlotlyJS.relayout!(fig, 
            scene2_xaxis_range = trace_point_layout[:scene_xaxis_range],
            scene2_yaxis_range = trace_point_layout[:scene_yaxis_range],
            scene2_zaxis_range = trace_point_layout[:scene_zaxis_range],
        )
    end
    PlotlyJS.add_trace!(fig, trace_line, row=1,col=1)
    PlotlyJS.add_trace!(fig, trace_3d, row=2,col=1)
    display(fig)
end

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
function plot_coordinate(vec_coord::Vector{Matrix{Float64}}, vec_point::Vector{Vector{Float64}})
    u,v,w =[],[],[]
    x,y,z = [], [], []
    # copy point 3 times
    foreach(1:3) do _
        x_tmp,y_tmp,z_tmp =
        map(1:3) do axis
            map(x->x[axis], vec_point)
        end   
        x = [x;x_tmp]
        y = [y;y_tmp]
        z = [z;z_tmp]
    end
    # the i is the basis x,y,z
    for i = 1:3
        u_tmp,v_tmp,w_tmp = 
        map(1:3) do axis
            map(M->M[axis,i], vec_coord)
        end
        u = [u;u_tmp]
        v = [v;v_tmp]
        w = [w;w_tmp]
    end
    PlotlyJS.plot(
        PlotlyJS.cone(
            x = x,
            y = y,
            z = z,
            u = u,
            v = v,
            w = w,
            sizemode = "absolute",
            sizeref = 1
        ),
        PlotlyJS.Layout(
            scene=PlotlyJS.attr(
                aspectmode="data",
            ),
            margin=PlotlyJS.attr(t=0, b=0, l=0, r=0)
        )
    )
end

function plot_coordinate(point::Vector{anten_point})
    local_coord = map(x->x.local_coord, point)
    position = map(i->[i.p...], point)
    plot_coordinate(local_coord, position)
end

export
    plot_point,
    plot_pattern,
    plot_E_vec_field,
    plot_E_vec_field_test,
    plot_pattern_2D,
    plot_coordinate,
    plot_pattern_all