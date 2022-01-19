module antenna

using LinearAlgebra
using Interpolations
using Chain
# pattern file
using CSV
using DataFrames
using Match
using Reexport

include("type.jl")
include("utils.jl")
include("plotting.jl")
@reexport using .type
@reexport using .plotting
@reexport using .utils

export
    cal_pattern,
    anten_read,
    rotate_vec_in_sph,
    rotate_vec_in_cart,
    rotate_pattern

# array function
# array factor
AF(point, θ, ϕ, k) = exp(1im * k * (point.x * sin(θ)cos(ϕ) + point.y * sin(θ)sin(ϕ) + point.z * cos(θ)))
# current source
Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k) ^ -1
# Iₛ(point, θₜ, ϕₜ, k) = 1
# pattern for one point
Pi(point, Pe, θ, ϕ, θₜ, ϕₜ, k) = Iₛ(point, θₜ, ϕₜ, k)Pe(θ, ϕ)AF(point, θ, ϕ, k)
# apply rotation
# XXX untest
rotate_vec_in_sph = (θ, ϕ, M) -> begin
    x, y, z = sph2cart(θ, ϕ, 1)
    vec = M * [x, y, z]
    θ1, ϕ1, .. = cart2sph(vec[1], vec[2], vec[3])
    # XXX not a elegen way to aviod lossing information
    # what if there is a rotation, boudary condition is when θ = 0
    # XXXX need a mapping function
    if θ == 0
        ϕ = ϕ
    elseif θ1 == 0
        ϕ = ϕ1
    end
    (θ1, ϕ1)
end
rotate_vec_in_cart = (vec,M) -> M*vec
vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0]


rotate_pattern = (pattern::anten_pattern, coord::Matrix{Float64}) -> begin

    θ_grid = [θ for θ in θ_default, ϕ in ϕ_default]
    ϕ_grid = [ϕ for θ in θ_default, ϕ in ϕ_default]

    # θ_grid_tmp = [rotate_vec_in_sph(θ,ϕ, coord')[1] for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    # ϕ_grid_tmp = [rotate_vec_in_sph(θ,ϕ, coord')[2] for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    # θ_grid_tmp = [θ for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    # ϕ_grid_tmp = [ϕ for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    # set_grid(θ_grid_tmp, ϕ_grid_tmp)
        
    Gθ_grid = [ begin
        θ,ϕ = rotate_vec_in_sph(θ,ϕ,coord');
        (coord*vec_θ(θ,ϕ))*pattern.θ(θ,ϕ) end
    for (θ,ϕ) = zip(θ_grid, ϕ_grid)]

    Gϕ_grid = [begin
        θ,ϕ = rotate_vec_in_sph(θ, ϕ, coord');
        (coord*vec_ϕ(θ,ϕ))*pattern.ϕ(θ,ϕ) end
    for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    
    vec_θ₁_map_grid = [vec_θ(θ,ϕ) for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    vec_ϕ₁_map_grid = [vec_ϕ(θ,ϕ) for (θ,ϕ) = zip(θ_grid, ϕ_grid)]

    Gθ = (dot.(Gθ_grid, vec_θ₁_map_grid))+
        (dot.(Gϕ_grid, vec_θ₁_map_grid))
    Gϕ = (dot.(Gϕ_grid, vec_ϕ₁_map_grid))+
        (dot.(Gθ_grid, vec_ϕ₁_map_grid))

    result = anten_pattern(
        θ = LinearInterpolation((θ_default, ϕ_default), Gθ),
        ϕ = LinearInterpolation((θ_default, ϕ_default), Gϕ)
    )
end

# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(point::Vector{anten_point}, θₜ, ϕₜ, k = k, θ = θ_default, ϕ = ϕ_default)

    # apply_rotation to the pattern
    # map((point::anten_point) ->
    #         point.pattern = rotate_pattern(point.pattern, point.local_coord),
    #     point
    # )

    # calculate result
    result_θ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1), size(point, 1))
    result_ϕ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1), size(point, 1))
    # Threads.@threads for (ind_point, point) = collect(enumerate(point))
    for (ind_point, point_i) = collect(enumerate(point))
        for (i, θ) = enumerate(θ), (j, ϕ) = enumerate(ϕ)
            result_θ[i, j, ind_point] = Pi(point_i.p, point_i.pattern.θ, θ, ϕ, θₜ, ϕₜ, k)
            result_ϕ[i, j, ind_point] = Pi(point_i.p, point_i.pattern.ϕ, θ, ϕ, θₜ, ϕₜ, k)
        end
    end

    result_θ = sum(result_θ, dims = 3)[:, :, 1] .|> abs
    result_ϕ = sum(result_ϕ, dims = 3)[:, :, 1] .|> abs

    anten_pattern(
        θ = LinearInterpolation((θ_default, ϕ_default), result_θ),
        ϕ = LinearInterpolation((θ_default, ϕ_default), result_ϕ),
    )
end





# order: the column number for θ,ϕ,Sθ,Sϕ, for HFSS is [2, 1, 4, 3]
function anten_read(filepath, type = "hfss"; unit = "db", order = [2, 1, 4, 3])
    #pick data from the dataframe
    unit_corrector = unit -> begin
        if unit == "db"
            x -> 10^(x / 10) |>sqrt
        else
            x -> x
        end
    end
    convert_from_hfss = x -> begin
        (
            θ = deg2rad.(x[:, order[1]]),
            ϕ = deg2rad.(x[:, order[2]]),
            Gθ = x[:, order[3]] .|> unit_corrector(unit),
            Gϕ = x[:, order[4]] .|> unit_corrector(unit)
        )
    end

    # read file
    raw_anten = CSV.File(filepath) |> DataFrame |>
                @match type begin
                    "hfss" => convert_from_hfss
                    _ => convert_from_hfss
                end

    _3Dvec_2_1Dvec = x -> LinRange(minimum(x), maximum(x), size(unique(x), 1))

    translate2Interpolation = (θ, ϕ, G) -> @chain begin
        reshape(G, (size(ϕ, 1), size(θ, 1)))
        permutedims((order[1], order[2]))
        LinearInterpolation((θ, ϕ), _)
    end
    θ = _3Dvec_2_1Dvec(raw_anten.θ)
    ϕ = _3Dvec_2_1Dvec(raw_anten.ϕ)

    anten_pattern(
        θ = translate2Interpolation(θ, ϕ, raw_anten.Gθ),
        ϕ = translate2Interpolation(θ, ϕ, raw_anten.Gϕ)
    )
end


# calculate result

# set_param(f=5.8e9)
# θₜ, ϕₜ = (0, 0) .|> deg2rad
# point = point_linear(N = 8, dx = 25e-3)
# set_point_loc_coord!.(point, [[1.0 0 0;0 1 0; 1 0 1]])
# P_norm = cal_pattern(point, θₜ, ϕₜ)  .|> scale2log |> normlog
# plot(P_norm[:,1])





# calculate interp
# pattern = anten_read("pattern.csv")
# plot(pattern, ϕ = pi / 2)




# apply rotation
# set_param(f=5.8e9)
# elem_pattern = anten_read("pattern.csv")
# point = point_linear(N = 8, dx = 25e-3, pattern = elem_pattern)
# set_point_loc_coord!.(point, [Matrix(1.0I, (3, 3))])
# #XXX FIXME wrong result
# θₜ, ϕₜ = (0, 90) .|> deg2rad
# global_pattern = @time cal_pattern(point, θₜ, ϕₜ)
# plot(global_pattern)
# # plot([global_pattern.θ(θ,ϕ) for θ=θ_default,ϕ=ϕ_default][:,271])


end
