module plotting

using ..type
using ..utils
using PlotlyJS
using LinearAlgebra
# θ_default, ϕ_default = (1:1:179, -180:1:179) .|> x -> deg2rad.(x)

function ploting_point(p::Vector{anten_point})
    expand_point = (vec_point) -> begin
        ([i.p.x for i = vec_point],
            [i.p.y for i = vec_point],
            [i.p.z for i = vec_point])
    end
    x, y, z = expand_point(p)
    plot(scatter(x = x, y = y, z = z,
        type = "scatter3d",
        marker = attr(
            size = 2
        ),
        mode = "markers"))
end
function ploting_point(p::Vector{Vector{Float64}})
    expand_point = (vec_point) -> begin
        ([i[1] for i = vec_point],
            [i[2] for i = vec_point],
            [i[3] for i = vec_point])
    end
    x, y, z = expand_point(p)
    plot(scatter(x = x, y = y, z = z,
        type = "scatter3d",
        marker = attr(
            size = 2
        ),
        mode = "markers"))
end
# ploting_point(point)

function ploting_pattern(pattern::anten_pattern)
    r = [pattern.θ(θ, ϕ)^2 + pattern.ϕ(θ, ϕ)^2 for (θ,ϕ) = zip(θ_grid,ϕ_grid)]
    θ = [θ for (θ,ϕ) = zip(θ_grid,ϕ_grid)]
    ϕ = [ϕ for (θ,ϕ) = zip(θ_grid,ϕ_grid)]
    tulple_matrix = sph2cart.(θ, ϕ, r)
    x = [i.x for i = tulple_matrix]
    y = [i.y for i = tulple_matrix]
    z = [i.z for i = tulple_matrix]

    plot(surface(x = x, y = y, z = z))

end
# ploting_pattern(global_pattern)

function plot_E_vec_field(pattern::anten_pattern; mode = "1vec")
    # get point from default θ and ϕ
    θ = [θ for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    ϕ = [ϕ for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    point_tuple_matrix = sph2cart.(θ, ϕ, 1)
    x = [i.x for i = vec(point_tuple_matrix)]
    y = [i.y for i = vec(point_tuple_matrix)]
    z = [i.z for i = vec(point_tuple_matrix)]

    vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
    vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0]
    # get pattern  E vector

    #convert E vector to 1D vector

    if mode == "2vec"
        vec_vec_θ = [vec_θ(θ, ϕ) * pattern.θ(θ, ϕ) for (θ,ϕ) = vec(zip(θ_grid, ϕ_grid))]
        vec_vec_ϕ = [vec_ϕ(θ, ϕ) * pattern.ϕ(θ, ϕ) for (θ,ϕ) = vec(zip(θ_grid, ϕ_grid))]
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
            vec_θ(θ, ϕ) * pattern.θ(θ, ϕ) .+  vec_ϕ(θ, ϕ) * pattern.ϕ(θ, ϕ)
            for (θ,ϕ) = zip(θ_grid, ϕ_grid)
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
    plot(
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

    # plot(scatter(x=x,y=y,z=z, 
    #             type="scatter3d",
    #              marker=attr(
    #                  size=2
    #             ),
    #              mode="markers"))
end

function plot_E_vec_field_test(vecE)
    # get point from default θ and ϕ
    θ = [θ for (θ,ϕ) = zip(θ_grid ,ϕ_grid)]
    ϕ = [ϕ for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    point_tuple_matrix = sph2cart.(θ, ϕ, 1)
    x = [i.x for i = vec(point_tuple_matrix)]
    y = [i.y for i = vec(point_tuple_matrix)]
    z = [i.z for i = vec(point_tuple_matrix)]

    vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
    vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0]
    # get pattern  E vector
    #convert E vector to 1D vector

        # vec_vec_globe = [[1,1,1] for (θ,ϕ) = zip(vec(θ_grid), vec(ϕ_grid))]
        u = [i[1] for i = vec(vecE)]
        v = [i[2] for i = vec(vecE)]
        w = [i[3] for i = vec(vecE)]
    plot(
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

    # plot(scatter(x=x,y=y,z=z, 
    #             type="scatter3d",
    #              marker=attr(
    #                  size=2
    #             ),
    #              mode="markers"))
end
export
    ploting_point,
    ploting_pattern,
    plot_E_vec_field,
    plot_E_vec_field_test
end