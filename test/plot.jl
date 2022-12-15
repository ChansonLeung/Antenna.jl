using LinearAlgebra
using Polymake
using Antenna
using Plots
using Lazy

plotlyjs()

# A = polytope.icosahedron()
# v = A.VERTICES |>eachrow .|>x->Float64.(collect(x[2:end])) 


# scatter(x=map(x->x[1],v), y=map(x->x[2],v), z=map(x->x[3],v), 
#         mode="markers", 
#         marker=attr(
#             size=4,
#             color="B",
#         ),
#         type="scatter3d") |> plot


dipole = read_hfss_pattern("test/dataset_accurate/rE_theta_phi_dipole.csv")
patch = read_hfss_pattern("test/elem/rE.csv")


@recipe function f(pattern::anten_pattern;seriestype=:heatmap)
    uv_grid =[(u,v) for u in LinRange(-1,1,300), v in LinRange(-1,1,300)] 
    u = @> uv_grid getindex.(1)
    v = @> uv_grid getindex.(2)
    θ = @. acos(1/sqrt(u^2+v^2+1))
    ϕ = @. atan(v,u)

    z = directivity(pattern).(θ, ϕ)  .|>x->10log10(x)

    seriestype --> seriestype
    z
end

@recipe function f(points::Vector{Float64};)
    seriestype-->:scatter3d
    (getindex.(points,1),
    getindex.(points,2),
    getindex.(points,3))
end


@recipe f(::Type{anten_point}, point::anten_point) = point.pattern
@recipe f(::Type{Vector{anten_point}}, points::Vector{anten_point}) = getfield.(points, :p)

points = point_rectangle(Nx = 10,Ny=10, dx=λ/2, dy=λ/2,pattern=dipole)
getfield.(points, :p)
points = getfield.(points, :p)

scatter3d(
    getindex.(points,1),
    getindex.(points,2),
    getindex.(points,3),
    markersize=2
)




# uv_grid_sort = sort(uv_grid, dims=1,lt=(l,r)->(l[1]<r[1]))
# PlotlyJS.plot(PlotlyJS.heatmap(z = getindex.(uv_grid_sort,2)))
