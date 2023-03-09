using PyCall
using PyFormattedStrings
using Antenna
import PlotlyJS
using GLMakie
using Tullio
using LoopVectorization

## ---------read pattern from hfss and get anten_pattern of patch---------
set_param(freq=2.7e9)
patch = read_hfss_pattern("test/elem/rE.csv")

## ---------create array for cylinder Vector{anten_point}---------
array_cylinder = create_array_cylinder(pattern = patch)

## ---------calculate result for cylinder---------
res_pattern, D = cal_pattern(array_cylinder, pi/2, 0, spin = true, optimal_directivity=:all)
plot_pattern(res_pattern)
point_cone(array_cylinder) |>PlotlyJS.plot

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

const u_const_grid = [u for u in range(-1,1,361), _ in range(-1,1,361)]
const v_const_grid = [v for _ in range(-1,1,361), v in range(-1,1,361)]
const θ_const_grid = (@. acos(1/(sqrt(u_const_grid^2+v_const_grid^2+1))))
const ϕ_const_grid = (@. atan(u_const_grid,v_const_grid))
function directivity_uv(
        pattern,
        θ = θ_const_grid,
        ϕ = ϕ_const_grid;
        power=100.0
    )
    directivity_beta(pattern, θ, ϕ, power=power) 
end

function sll_gain(pattern; power=100.0, min_distance=4)
    #FIXME how to convert fastly
    dir = directivity_uv(pattern, power = power);
    idx_maxima = py"peak_local_max"(dir, min_distance=min_distance)
    SLL = map(x-> dir[x[1],x[2]], eachrow(idx_maxima))
    sort!(SLL, order=Base.Order.Reverse)
    max_sll = 10log10(SLL[2]/SLL[1])
    maxgain= 10log10(SLL[1] )
    max_sll,maxgain
end

iter_data = [@SVector[0.0,0.0,0.0]]
best_solution = []
best_fit = [0.0]
f_count = 0
function objective(x, array=array_plane)
    global f_count += 1
    ## --------- calculate pattern with input x---------
    for idx in eachindex(array)
        array[idx].coeffi = x[idx]
    end
    res_pattern, D = cal_pattern(array, 0, 0)
    
    ##--------- calculate max SLL and gain by using scikit-image's peak_local_max ---------
    max_sll, maxgain = sll_gain(res_pattern, power=sum(abs.(D)))
    # dir = directivity_uv(res_pattern, power = sum(abs.(D)));
    # idx_maxima = py"peak_local_max"(dir, min_distance=4)
    # SLL = map(x-> dir[x[1],x[2]], eachrow(idx_maxima))
    # sort!(SLL, order=Base.Order.Reverse)
    # max_sll = 10log10(SLL[2]/SLL[1])
    # maxgain= 10log10(SLL[1] )

    ## ---------save important data ---------
    push!(iter_data, @SVector[Float64(f_count) , max_sll, maxgain])
    if max_sll < best_fit[end]
        push!(best_fit, max_sll)
        push!(best_solution, x)
    end

    ## ---------update plot and print result---------
    sc1[1] = res_pattern
    sc2[2] = abs.(D);
    println(f"{ f_count:.1f }: { max_sll:.2f }dB,   { maxgain:.2f }dB")
    max_sll
end

## ---------profile objective function---------
input = rand(400);
@time objective(input, array_plane);
@profview objective(input, array_plane);

## --------- optimizing! ---------
using BlackBoxOptim
bboptim_res =  bboptimize(objective;SearchRange=(0,1),TargetFitness=-30.0, NumDimensions=400, MaxSteps=500_000)

## --------- plot fitness line---------
fig = Figure(viewmode=:data)
ax = Axis(fig[1,1])
tt = (getindex.(iter_data,2)) |> x->replace(x, -Inf=>-15)
lines!(ax, getindex.(iter_data,2))
fig

## ---------plot local max in heatmap
x_end = best_solution[end]
for (p_i, xi) in zip(array_plane, x_end)
    p_i.coeffi = xi
end
res, D = cal_pattern(array_plane, 0, 0);
dir = directivity_uv(res).|> x->10log10(x)
fig = heatmap(dir)
scatter!((py"peak_local_max"(dir, min_distance= 2)), color=:pink)
idx_maxima = py"peak_local_max"(dir, min_distance= 2);
fig

##saveing
# using JLD2
# iter_data
# best_solution
# best_fit
# jldsave("topic/20x20plane_lowsidelope.jld2";array_plane, iter_data, best_solution, best_fit)
# ttt = load("topic/20x20plane_lowsidelope.jld2")