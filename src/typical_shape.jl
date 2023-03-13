using Antenna
using LinearAlgebra
using StaticArrays
using Rotations
export 
    create_array_cylinder,
    create_array_sphere,
    point_cone,
    # mesh
    create_array_from_meshes_sample,
    create_facets_from_file,
    rotate_facet,
    vecfacet2ply
    
## ---------define array create function---------
function create_array_cylinder(; R=2λ, dh=λ / 2, N_layers=nothing, N_in_eachlayer=nothing, pattern=pattern_identity)
    N_in_eachlayer === nothing && (N_element_in_eachlayer = ceil(Int, 2pi * R / (λ / 2)))
    N_layers === nothing && (N_layers = ceil(Int, 2R / (λ / 2)))
    map(range(0, 2pi - 2pi * (1 / N_element_in_eachlayer), N_element_in_eachlayer)) do ϕ
        map(range(0, N_layers * dh, N_layers)) do h
            x = R * cos(ϕ)
            y = R * sin(ϕ)
            z = h
            local_coord = zeros(3, 3)
            u, v, w = eachcol(local_coord)
            u .= [0, 0, -1]
            w .= [cos(ϕ), sin(ϕ), 0]
            v .= cross(w, u)
            anten_point(p=[x, y, z], local_coord=local_coord, pattern=pattern)
        end
    end |> x -> vcat(x...)
end

function point_cone(array)
    x = getfield.(array, :p) .|> x->x[1]
    y = getfield.(array, :p) .|> x->x[2]
    z = getfield.(array, :p) .|> x->x[3]
    coords = getfield.(array, :local_coord)
    u = map(x->x[1,1], coords)
    v = map(x->x[2,1], coords)
    w = map(x->x[3,1], coords)
    PlotlyJS.cone(
        x = x, y = y, z = z,
        u = u, v = v, w = w,
        showscale=false,
        colorscale=[[0,PlotlyJS.colors.Oranges[5]], [1,PlotlyJS.colors.Oranges[5]]],
    )
end

## --------- define circle create---------

## create sphere
function create_array_sphere(
    arc = λ / 2*1.2,
    R = 3λ;
    pattern = pattern_identity
)

    #plane-plane intersect, find the right vector of element, and calculate the RotZ  angle
    function RotZ_rot_angle(ϕ, elem_vec_z, elem_vec_x)
        vec_polorazing = [-sin(ϕ), cos(ϕ), 0]
        elem_vec_x_aims = cross(vec_polorazing,elem_vec_z) /norm(elem_vec_z)
        # acos(clamp(dot(elem_vec_x, elem_vec_x_aims), -1,1))
        -sign(dot(cross(elem_vec_x, elem_vec_z), elem_vec_x_aims)) * atan(norm(cross(elem_vec_x,elem_vec_x_aims)),dot(elem_vec_x, elem_vec_x_aims))
    end

    map(LinRange(0, pi, ceil(Int64, pi*R/(arc)))) do θ
        r = sin(θ) * R
        div = ceil(Int64, 2pi * r / (arc))
        num_per_ring = ceil(Int64, 2pi * r / (arc))

        # div > 1 && map(LinRange(0.0, 2pi*(num_per_ring - 1)/num_per_ring , num_per_ring-1)) do ϕ
        ringϕ = div >1 ? LinRange(0.0, 2pi*(num_per_ring - 1)/num_per_ring , num_per_ring-1) : [0.]
        map(ringϕ) do ϕ
            local_coord = zeros(3,3)
            u,v,w =  map(i->(@view local_coord[:, i]),1:3)
            u .= vec_θ(θ, ϕ) 
            v .= vec_ϕ(θ, ϕ)  
            w .= cross(u,v)  
            # rotate to the same polarization, with polarization vector is phi = phi
            local_coord *= RotZ(RotZ_rot_angle(deg2rad(0), w, u))
            p = anten_point(p=[r * cos(ϕ), r * sin(ϕ), R*cos(θ)], pattern=pattern, local_coord=local_coord)
        end
    end |>  x->reduce(vcat,x)|> x->filter(i->isa(i ,anten_point) ,x) |>Vector{anten_point}
end


##--------- conformal shpae generate ---------
struct myfacet
    normal::Vector{Float64}
    vertex::Vector{Vector{Float64}}
end
xyzString2Matrix(str) = str |> String |> IOBuffer |> x -> CSV.read(x, DataFrame) |> Matrix
Matrix2xyzString(mat) = replace(mat |> string, "[" => "", "]" => "", ";" => "\n")
#plane-plane intersect, find the right vector of element, and calculate the RotZ  angle
function RotZ_rot_angle(ϕ, elem_vec_z, elem_vec_x)
    vec_polorazing = [-sin(ϕ), cos(ϕ), 0]
    elem_vec_x_aims = cross(vec_polorazing,elem_vec_z) /norm(elem_vec_z)
    # acos(clamp(dot(elem_vec_x, elem_vec_x_aims), -1,1))
    -sign(dot(cross(elem_vec_x, elem_vec_z), elem_vec_x_aims)) * atan(norm(cross(elem_vec_x,elem_vec_x_aims)),dot(elem_vec_x, elem_vec_x_aims))
end
function find_facet(mesh, d, point)
    s = KBallSearch(mesh, 1, MetricBall(d))
    idx_points = search(point, s)
    map(idx_points) do i
        mesh.topology.connec[i]
        # mesh.vertices[i]
    end
end
function create_facets_from_file(file_name="test/dataset/surface.STL")
    vec_facet = Vector{myfacet}()
    println("-----------")
    open(file_name, "r") do f
        while true
            line = readline(f)
            if line == ""
                break
            end
            facet_i = myfacet(zeros(3), fill([0, 0, 0], 3))
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
                facet_i.vertex[i][1] = parse(Float64, match_res[1])/1000
                facet_i.vertex[i][2] = parse(Float64, match_res[2])/1000
                facet_i.vertex[i][3] = parse(Float64, match_res[3])/1000
            end
            readline(f) # endloop
            readline(f) # end facet
            push!(vec_facet, facet_i)
        end
    end
    vec_facet
end

function rotate_facet(vec_facet)
    for facet_i in (vec_facet)
        for v_i in Ref.(facet_i.vertex)
            v_i[] .= RotX(pi / 2)*v_i[]
        end
    end
end

function vecfacet2ply(vec_facet)
    p = map(x -> x.vertex, vec_facet) |> x -> reduce(vcat, x) |> x -> reduce(vcat, x')'
    vertex = eachcol(p) |>collect 
    dict_vertex = map((i, idx) -> i => idx, vertex, eachindex(vertex)) |> Dict{Vector{Float64}, Int64}
    facet_connections = map(getfield.(vec_facet, :vertex)) do (v1,v2,v3)
        connect((dict_vertex[v1], dict_vertex[v2], dict_vertex[v3]))
    end
    plane_mesh = SimpleMesh(vertex .|> Meshes.Point3, facet_connections)
end



function create_array_from_meshes_sample(point_sample, mesh; pattern = pattern_identity)
    ## create array
    plane_vec_norm = map(point_sample) do p
        mesh = mesh
        points_idx = find_facet(mesh, λ/2, p)[1].indices |> x->[x...] 
        points = points_idx .|> x->mesh.vertices[x] .|> x->[x.coords.coords...]
        vec_norm = cross(points[2]-points[1], points[3]-points[1])
        # FIXME unknow reason for the minus
        vec_norm = -vec_norm./norm(vec_norm)
    end
    plane_array = map(point_sample, plane_vec_norm) do p, vec_norm
        θ,ϕ,r = cart2sph(p.coords...)
        local_coord = zeros(3,3)
        u,v,w =  map(i->(@view local_coord[:, i]),1:3)
        w .= vec_norm
        v .= cross([0,0,1.],w) |>x->x./norm(x)
        u .= cross(v, w)
        local_coord *= RotZ(RotZ_rot_angle(deg2rad(-90), w, u))
        anten_point(p=p.coords, pattern=pattern, local_coord=local_coord)
    end
end

