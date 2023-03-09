using Antenna
using LinearAlgebra
using StaticArrays
using Rotations
export 
    create_array_cylinder,
    create_array_sphere,
    point_cone
    
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
    R = 3λ
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
            p = anten_point(p=[r * cos(ϕ), r * sin(ϕ), R*cos(θ)], pattern=patch, local_coord=local_coord)
        end
    end |>  x->reduce(vcat,x)|> x->filter(i->isa(i ,anten_point) ,x) |>Vector{anten_point}
end