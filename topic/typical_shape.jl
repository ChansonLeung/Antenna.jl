using Antenna
using GLMakie
using LinearAlgebra
using PyFormattedStrings
using StaticArrays
using PyCall

## ---------define array create function---------
function create_array_cylinder(; R=5λ, dh=λ / 2, N_layers=nothing, N_in_eachlayer=nothing, pattern=pattern_identity)
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

## ---------visualize with GLMakie---------
include("./makie_antenna.jl")


# surface!(x,y,z,color=r, colormap="lightrainbow")
# create_figure()

## ---------read pattern from hfss and get anten_pattern of patch---------
set_param(freq=2.7e9)
patch = read_hfss_pattern("test/elem/rE.csv")
## ---------create array for cylinder Vector{anten_point}---------
array_cylinder = create_array_cylinder(pattern = patch)

## ---------calculate result for cylinder---------
res_pattern, D = cal_pattern(array_cylinder, pi/2, 0, spin = true, optimal_directivity=:all)
plot_pattern(res_pattern)
plot_point(array_cylinder)
# point_cone(array_cylinder) |>PlotlyJS.plot

## ---------create plane array---------
array_plane = point_rectangle(Nx = 20, Ny=20, dx=λ/2, dy=λ/2, pattern=patch)
## ---------calculate result for plane---------
res_pattern, D = cal_pattern(array_plane, 0,0)
plot_pattern(res_pattern)
## --------- create array for sphere with ring style

##--------visualize in GLMakie---------
fig = Figure(resolution=(1200,600), viewmode=:fit, aspect=:data) 
display(fig)
ax1, sc1 = pattern3d(fig[1,2], res_pattern);
ax2, sc2 = arraypoints(fig[1,1], array_plane, zeros(size(array_plane)))
# slider_θ = Slider(fig[2,:], range=-180:180, startvalue=0)
# slider_ϕ = Slider(fig[3,:], range=-180:180, startvalue=0)
# lift(slider_θ.value, slider_ϕ.value) do θ,ϕ
#     θ,ϕ = mod_angle_rad(deg2rad(θ),deg2rad(ϕ))
#     @time res_pattern, D = cal_pattern(array_plane, θ,ϕ, spin = true, optimal_directivity=:all);
#     sc1[1] = res_pattern
#     sc2[2] = abs.(D);
#     sc3[2] = angle.(D);
# end
# fig
##---------Optimization---------
using Optim
using PyCall
using StaticArrays

py"""
from skimage.feature import peak_local_max
"""
final_result = nothing
exit_flag = false

final_result = nothing
function advanced_time_control(x)
    global final_result = x
    exit_flag
end 
Threads.@spawn begin
    read(stdin, 1)
    println("end optimization")
    global exit_flag = true
end

function directivity_uv(pattern, 
        u = [u for u in range(-1,1,181), _ in range(-1,1,181)],
        v = [v for _ in range(-1,1,181), v in range(-1,1,181)],
        θ = (@. acos(1/(sqrt(u^2+v^2+1)))),
        ϕ = (@. atan(u,v)),
    )
    # directivity_beta(pattern, θ, ϕ) 

    @. 4pi*1 / 2(120pi) * (abs(pattern.θ(θ, ϕ))^2 + abs(pattern.ϕ(θ, ϕ))^2)/100.0
end
function sll_gain(pattern)
    #FIXME how to convert fastly
    dir = directivity_uv(pattern);
    idx_maxima = py"peak_local_max"(dir, min_distance=3)
    SLL = map(x-> dir[x...],eachrow(idx_maxima))
    sort!(SLL, order=Base.Order.Reverse)
    SLL[2]/SLL[1], SLL[1] 
end

iter_data = [@SVector[0.0,0.0,0.0]]
best_solution = []
best_fit = [0.0]
global f_count = 0
function f(x)
    global f_count += 1
    for idx in eachindex(array_plane)
        array_plane[idx].coeffi = x[idx]
    end
    res_pattern, D = cal_pattern(array_plane, 0, 0)
    
    # sc1[1] = res_pattern
    # sc2[2] = abs.(D);
    # ---------
    # max_sll, maxgain = sll_gain(res_pattern)
    dir = directivity_uv(res_pattern);
    idx_maxima = py"peak_local_max"(dir, min_distance=3)
    SLL = map(x-> dir[x...],eachrow(idx_maxima))
    sort!(SLL, order=Base.Order.Reverse)
    max_sll = SLL[2]/SLL[1]
    maxgain= SLL[1] 
    #--------- 


    # push!(iter_data, @SVector[Float64(f_count) , max_sll, maxgain])

    ssl = 10log10(max_sll)
    maxgain = 10log10(maxgain)

    if ssl < best_fit[end]
        push!(best_fit, max_sll)
        push!(best_solution, x)
    end
    println(f"{ f_count:.1f }: { ssl:.2f }dB,   { maxgain:.2f }dB")
    ssl
end
# input=ones(256)
# res = optimize(f, input, ParticleSwarm(lower=zeros(256), upper=ones(256)), Optim.Options(callback = advanced_time_control))
# res_pattern, D = cal_pattern(array_plane, 0,0);
# sc1[1] = res_pattern
# sc2[2] = abs.(D);
# SLL_maxgain(res_pattern)[1] |> x->10log10(x)
## --------- bloack box---------
using BlackBoxOptim

bboptimize(f;SearchRange=(0,1),TargetFitness=-30.0, NumDimensions=400, MaxSteps=500_000)


## ---------test local maxima---------
input = rand(400);
@profview f(input);
@time res_pattern, D = cal_pattern(array_plane, 0, 0);
@benchmark sll_gain($res_pattern)

fig = Figure(viewmode=:data)
ax = Axis(fig[1,1])
tt = (getindex.(iter_data,2).|> x->10log10(x)) |> x->replace(x, -Inf=>-15)
lines(tt)
lines!(ax, getindex.(iter_data,2).|> x->10log10(x))
iter_data


# using JLD2
# jldsave("" )
# load("")

# final_result.value
x_end = best_solution[end]
best_solution
best_fit
for (p_i, xi) in zip(array_plane, x_end)
    p_i.coeffi = xi
end
@time res, D = cal_pattern(array_plane, 0, 0);
plot_pattern(res)
dir = directivity_uv(res).|> x->10log10(x)
heatmap(dir)
scatter!((py"peak_local_max"(dir, min_distance= 3)))
@time idx_maxima = py"peak_local_max"(dir, min_distance= 3);

SLL = map(x-> dir[x...],eachrow(idx_maxima))
PlotlyJS.heatmap(z=directivity_uv(res).|> x->10log10(x)) |> PlotlyJS.plot

@time f(x_end);
@time sll_gain(res);

10log10(t1)