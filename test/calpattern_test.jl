##
using Antenna
using BenchmarkTools
using Tullio
using LinearAlgebra
using PlotlyJS

patch = read_hfss_pattern("test/elem/rE.csv")
set_param(freq=2.7e9)

## create sphere
arc = λ / 2
R = 3λ
#plane-plane intersect, find the right vector of element, and calculate the RotZ  angle
function RotZ_rot_angle(ϕ, elem_vec_z, elem_vec_x)
    vec_polorazing = [-sin(ϕ), cos(ϕ), 0]
    elem_vec_x_aims = cross(vec_polorazing,elem_vec_z) /norm(elem_vec_z)
    # acos(clamp(dot(elem_vec_x, elem_vec_x_aims), -1,1))
    -sign(dot(cross(elem_vec_x, elem_vec_z), elem_vec_x_aims)) * atan(norm(cross(elem_vec_x,elem_vec_x_aims)),dot(elem_vec_x, elem_vec_x_aims))
end

array_sphere = map(LinRange(0, pi, ceil(Int64, pi*R/(arc)))) do θ
    r = sin(θ) * R
    div = ceil(Int64, 2pi * r / (arc))
    num_per_ring = ceil(Int64, 2pi * r / (arc))
    div > 1 && map(LinRange(0.0, 2pi*(num_per_ring - 1)/num_per_ring , num_per_ring-1)) do ϕ
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

###create vector field
using Rotations

x = getfield.(array_sphere, :p) .|> x->x[1]
y = getfield.(array_sphere, :p) .|> x->x[2]
z = getfield.(array_sphere, :p) .|> x->x[3]
# rot
urot,vrot,wrot = begin 
    local x_t,y_t,z_t =RotY(-pi/2)*[x y z]'  |>eachrow .|>collect
    local θ, ϕ, r=   cart2sph.(x_t,y_t,z_t) |> x-> (1:3 .|> i->getindex.(x,i))

    local u = @. cos(θ)cos(ϕ)
    local v = @. cos(θ)sin(ϕ)
    local w = @. -sin(θ)
    RotY(pi/2)*[u v w]'|>eachrow .|>collect
end
# norot
cos_vec(v1,v2) =acos(dot(v1,v2)/(norm(v1)norm(v2)))
u_nor, v_nor, w_nor = begin 
    local θ,ϕ, r=   cart2sph.(x,y,z) |> x-> (1:3 .|> i->getindex.(x,i))
    local u = @. cos(θ)cos(ϕ)
    local v = @. cos(θ)sin(ϕ)
    local w = @. -sin(θ)
    [u,v,w]
end
# rot by local_coord
u_coo, v_coo, w_coo = begin 
    coords = getfield.(array_sphere, :local_coord)
    local u = map(x->x[1,1], coords)
    local v = map(x->x[2,1], coords)
    local w = map(x->x[3,1], coords)
    [u, v, w]
end
# u,v,w = urot,vrot, wrot
# u,v,w = u_nor,v_nor, w_nor
u,v,w = u_coo,v_coo, w_coo

plot([cone(
    x = x,y=y,z=z,u=u,v=v,w=w,
    showscale=true,
    colorscale="Jet"
), trace_surface])

## calculate

res_pattern, D = @time cal_pattern(array_sphere, deg2rad(0.), 0.0, spin=true, optimal_directivity=:all);
# res_pattern, D = @time cal_pattern(array_sphere, deg2rad(0.), 0.0, spin=true);
plot_pattern(res_pattern)

# E_t = begin 
#     map(array_sphere) do i
#         (linear_interp_pattern_2com(getfield(i,:pattern_grid)).θ(0.,0.) ,
#         linear_interp_pattern_2com(getfield(i,:pattern_grid)).ϕ(0.,0.) ) 
#     end
# end
# filter(x->abs(angle(x[1]) - angle(x[2])) > deg2rad(20) && !(abs(x[1]) > abs(x[2])*100), E_t) .|>
# x->(angle(x[1]), angle(x[2]), abs(x[1]), abs(x[2]))

## plot sphere and point
# sphere grid
p = [0.975R.*[sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)] for θ in LinRange(0,pi,200) , ϕ in LinRange(-pi,pi,200)] 

trace_surface = begin 
    local x,y,z = [getindex.(p, i) for i in 1:3]
    surface(
        x = x, y = y, z = z,
        # surfacecolor = ones(Float64,length(z)),
        colorscale="Greys",
        surfacecolor = zeros(size(x)...),
        # lighting=attr(ambient=0.4, diffuse=0.5, roughness = 0.9, specular=0.6, fresnel=0.2),
        lightposition = attr(x=1000,y=1000,z=0),
        showscale=false,
    )
end
trace_point = begin
    local x,y,z = [map(anten_point->anten_point.p[i], array_sphere) for i in 1:3]
    scatter3d(
        x = vec(x), y = vec(y), z = vec(z), 
        mode="markers",
        marker = attr(size=10,color=(abs.(D).^2 |> x->x./maximum(x)),colorscale="Heat", showscale=true,
            colorbar=attr(
                tickfont=attr(
                    size=30
                )

            )
        ),
        # marker = attr(size=3,color=(abs.(D).^2 ),colorscale="Heat", showscale=true,)
    )
end



maximum_directivity(res_pattern)[2] .|> rad2deg

fig = plot(
    [trace_point, trace_surface],
    Layout(
        scene = attr(
            xaxis=attr(
                showgrid=false,
                showticklabels = false,
                showbackground = false,
                showtickprefix = true,
                title=""
            ),
            yaxis=attr(
                showgrid=false,
                showticklabels = false,
                showbackground = false,
                title=""
            ),
            zaxis=attr(
                showgrid=false,
                showticklabels = false,
                showbackground = false,
                title=""
            ),
            camera = attr(
                eye=attr(x=0.8, y=0.8, z=0.8)
            ),
        ),
        autosize=false,
        width=800, height=600,
    )
)
# display(fig)
# savefig(fig, "./microwaveweek/1.png",width=800, height=600, )
        

## plot quiver
