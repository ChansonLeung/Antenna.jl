using Revise
using Antenna
using BenchmarkTools
using Evolutionary
using StaticArrays

const target = (0.0, 0.0)
set_param(f=24e9)
const c = 299792458
const f = 24e9
const λ = c / f

function update_point!(array, dx̂, Î)
    N = length(array)
    for (index_l, index_c) in zip(LinearIndices(array), CartesianIndices(array))
        array[index_c].coeffi = Î[index_l] * exp(1im * 2pi * Î[index_l + N])
    end
end

function update_point_linear!(array, dx̂, Î)
    N = length(dx̂)
    dx̂[1] = dx̂[1]/2
    px = accumulate(+, dx̂)
    for index_l in LinearIndices(zeros(N))
        i = index_l
        I = 1* exp(1im * Î[i] * pi)
        array[index_l] = anten_point(p=[px[i], 0, 0], pattern=pattern_identity,     coeffi=I)
    end
end

function init_point(Nx,Ny)
    array = anten_point[]
    for i in 1:Nx, j in 1:Ny
        push!(array, anten_point(p=[λ/2*i, λ/2*j, 0], pattern=pattern_identity, coeffi=1))
    end
    array
end

count = 0
SLL_last_max = 0
SLL_mark = []

function transform(array, x;debug=false)
    global count = count +1

    update_point!(array, λ/2 * (array |>length|> ones) , x)
    res = cal_pattern(array,target...)

    if debug
        res
    else
        SLL,_ = SLL_maxgain(res) .|> x->10log10(x)
        if SLL > SLL_last_max
           global SLL_mark = x
           global SLL_last_max = SLL
           @show SLL, count 
        end
        (SLL - 25)^2
    end
end


array = init_point(10,10)

result = Evolutionary.optimize(
            x->transform(array, x),
            BoxConstraints(zeros(200),ones(200)),
            ones(200),
            GA(populationSize = 50, selection = susinv,
                crossover = DC, mutation = PLM()),
            Evolutionary.Options(parallelization=:thread, iterations = 10_000)
)

# p = transform(array, SLL_mark,debug=true)
# p = cal_pattern(array, 0.0,0.0, k=k,θ = θ_default, ϕ = LinRange(-180., 180., 361) .|> deg2rad)
# plot_pattern(p, θ = θ_default, ϕ = LinRange(-180., 180., 361).|> deg2rad)

# Evolutionary.minimizer(result) |> x-> transform(array, x,debug=true) |> plot_pattern_2D
# result
# @show Evolutionary.minimizer(result) 


# # speed test
# function test()
#     array = init_point(20*20)
#     for i in 1:100
#         update_point!(array, λ/2 * ones(20) , ones(20))
#         res =@time cal_pattern(array,target...)
#     end
# end
# @time test()

# test single

# array = init_point(20*20)
# update_point!(array, λ/2 * ones(20) , ones(20))
# res =@time cal_pattern(array,target...)
# plot_pattern_2D(res)