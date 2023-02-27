## convert stl to vec_facet
using CSV
using DataFrames
using Meshes
using Rotations

struct facet
    normal::Vector{Float64}
    vertex::Vector{Vector{Float64}}
end

xyzString2Matrix(str) = str |> String |> IOBuffer |> x -> CSV.read(x, DataFrame) |> Matrix
Matrix2xyzString(mat) = replace(mat |> string, "[" => "", "]" => "", ";" => "\n")

vec_facet = Vector{facet}()
println("-----------")
open("test/dataset/surface.STL", "r") do f
    while true
        line = readline(f)
        if line == ""
            break
        end
        facet_i = facet(zeros(3), fill([0, 0, 0], 3))
        match_res = match(r" *facet normal (.*) (.*) (.*)", line)

        if match_res === nothing
            line != "solid surface" && line != "endsolid" && throw("convert error line: " * line)
            continue
        else
            facet_i.normal[1] = parse(Float64, match_res[1])
            facet_i.normal[2] = parse(Float64, match_res[2])
            facet_i.normal[3] = parse(Float64, match_res[3])
        end
        readline(f) # outer loop
        for i in 1:3
            line = readline(f) # vertex
            match_res = match(r" *vertex (.*) (.*) (.*)", line)
            facet_i.vertex[i][1] = parse(Float64, match_res[1])
            facet_i.vertex[i][2] = parse(Float64, match_res[2])
            facet_i.vertex[i][3] = parse(Float64, match_res[3])
        end
        readline(f) # endloop
        readline(f) # end facet
        push!(vec_facet, facet_i)
    end
end

# rotate point for human readable
for facet_i in (vec_facet)
    # facet_i.vertex[1] = [0,0,0]
    for v_i in Ref.(facet_i.vertex)
        v_i[] .= RotX(pi / 2)*v_i[]
    end
end
vec_facet

## plot facets
p = map(x -> x.vertex, vec_facet) |> x -> reduce(vcat, x) |> x -> reduce(vcat, x')'

p_normal = map(x -> x.normal, vec_facet) |> x -> reduce(vcat, x')'
x_raw = p[1, :]
y_raw = p[2, :]
z_raw = p[3, :]
u = p_normal[1, :] |> x -> repeat(x, inner=3)
v = p_normal[2, :] |> x -> repeat(x, inner=3)
w = p_normal[3, :] |> x -> repeat(x, inner=3)

# mesh3d(x_raw, y_raw, z_raw, alpha=0.9, camera=(45, 45))
# scatter3d!(x_raw, y_raw, z_raw, ms=0.4)

## convert stl to ply
vertex = eachcol(p) |>collect
dict_vertex = map((i, idx) -> i => idx, vertex, eachindex(vertex)) |> Dict{Vector{Float64}, Int64}

facet_connections = map(getfield.(vec_facet, :vertex)) do (v1,v2,v3)
    connect((dict_vertex[v1], dict_vertex[v2], dict_vertex[v3]))
end

mesh = SimpleMesh(vertex .|> Meshes.Point3, facet_connections)
# fig, scene = viz(mesh)
# viz!(points_sample,color=:red, size=12)

## create array
import GLMakie
import MeshViz
# fig, scene = MeshVizviz(mesh)
MeshViz.viz(sphere_disc,show_facets=true)
sphere_disc = discretize(Sphere(Point3(0,0,0), R),RegularDiscretization(100,100) )
sphere_sample = sample(sphere_disc, MinDistanceSampling(λ/2*1.01, ρ=0.9)) |>collect

##
poisson_array_position =  sphere_sample.|> x->[x.coords.coords...] 
# poisson_array_position = filter(x->x[3]>0, poisson_array_position)
 
poisson_array = map(poisson_array_position) do p
    θ,ϕ,r = cart2sph(p...)
    local_coord = zeros(3,3)
    u,v,w =  map(i->(@view local_coord[:, i]),1:3)
    u .= vec_θ(θ, ϕ) 
    v .= vec_ϕ(θ, ϕ)  
    w .= cross(u,v)  
    local_coord *= RotZ(RotZ_rot_angle(deg2rad(0), w, u))
    anten_point(p=p, pattern=patch, local_coord=local_coord)
end

poi_u, poi_v, poi_w = begin 
    local coords = getfield.(poisson_array, :local_coord)
    local u = map(x->x[1,1], coords)
    local v = map(x->x[2,1], coords)
    local w = map(x->x[3,1], coords)
    [u,v,w]
end
poi_x,poi_y,poi_z = begin 
    local x = getindex.(poisson_array_position, 1)
    local y = getindex.(poisson_array_position, 2)
    local z = getindex.(poisson_array_position, 3)
    [x,y,z]
end
plot([cone(
    x = poi_x, y = poi_y, z = poi_z,
    u = poi_u, v = poi_v, w = poi_w,
    showscale=true,
    colorscale="Jet"
), trace_surface])
plot([scatter3d(
    x = poi_x, y = poi_y, z = poi_z,
    mode="markers",
    marker = attr(size=3,colorscale="Heat", showscale=true,
        colorbar=attr(
            tickfont=attr(
                size=30
            )

        )
    ),
    # u = poi_u, v = poi_v, w = poi_w,
    showscale=true,
    colorscale="Jet"
), trace_surface])
##
res_pattern, D = @time cal_pattern(filter(x->x.p[3] > 0, poisson_array), deg2rad(0.), 0.0, spin=true, optimal_directivity=:all);


plot_pattern(res_pattern)
plot_point(poisson_array)
max_dirs_poisson = map(-90:90 .|> deg2rad .|>x->mod_angle_rad(x,0.)[1], -90:90 .|> deg2rad .|>x->mod_angle_rad(x,0.)[2]) do θ,ϕ
    println(rad2deg(θ)," ", rad2deg(ϕ)) 
    res_pattern, D =  cal_pattern(poisson_array, θ, ϕ, spin=true, optimal_directivity=:all);
    maximum_directivity(res_pattern)
end


fig  = plot([
    scatter(x=-90:90, y = map(x->x[1],max_dirs_poisson) |>x->10log10.(x), name="poisson"), 
    scatter(x=-90:90, y = map(x->x[1],max_dir_ring) |>x->10log10.(x), name="ring")
])
# savefig(fig, "./microwaveweek/gain.html")

# using JLD2
# save("./microwaveweek/poisson_points.jld2", "data",reduce(hcat,poisson_array_position))
# load("./microwaveweek/poisson_points.jld2")["data"]

# poisson_array[1]



