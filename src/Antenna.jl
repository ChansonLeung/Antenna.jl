module Antenna

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

export
    cal_pattern,
    anten_read,
    rotate_vec_in_sph,
    rotate_vec_in_cart,
    rotate_pattern,
    directivity,
    radiation_intensity,
    radiated_power,
    gain,
    Iₛ

# array function
# array factor
AF(point, θ, ϕ, k) = exp(1im * k * (point.x * sin(θ)cos(ϕ) + point.y * sin(θ)sin(ϕ) + point.z * cos(θ)))
# current source
# Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k) ^ -1
Iₛ(point, θₜ, ϕₜ, k) = 1
# Iₛ(point, θₜ, ϕₜ, k) = 1
# pattern for one point
Pi(point, Pe, θ, ϕ, θₜ, ϕₜ, k) = Pe(θ, ϕ)AF(point, θ, ϕ, k)
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
    
    θ_grid::Matrix{Float64} = [θ for θ in θ_default, ϕ in ϕ_default]
    ϕ_grid::Matrix{Float64} = [ϕ for θ in θ_default, ϕ in ϕ_default]
        
    Gθ_grid::Matrix{Vector{ComplexF64}} = [ begin
        θ,ϕ = rotate_vec_in_sph(θ,ϕ,coord');
        (coord*vec_θ(θ,ϕ))*pattern.θ(θ,ϕ) end
    for (θ,ϕ) = zip(θ_grid, ϕ_grid)]

    Gϕ_grid::Matrix{Vector{ComplexF64}} = [begin
        θ,ϕ = rotate_vec_in_sph(θ, ϕ, coord');
        (coord*vec_ϕ(θ,ϕ))*pattern.ϕ(θ,ϕ) end
    for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    
    vec_θ₁_map_grid::Matrix{Vector{Float64}} = [vec_θ(θ,ϕ) for (θ,ϕ) = zip(θ_grid, ϕ_grid)]
    vec_ϕ₁_map_grid::Matrix{Vector{Float64}} = [vec_ϕ(θ,ϕ) for (θ,ϕ) = zip(θ_grid, ϕ_grid)]

    Gθ = (dot.(Gθ_grid, vec_θ₁_map_grid)) .+
        (dot.(Gϕ_grid, vec_θ₁_map_grid))
    Gϕ = (dot.(Gϕ_grid, vec_ϕ₁_map_grid)) .+
        (dot.(Gθ_grid, vec_ϕ₁_map_grid))

    result = anten_pattern(
        θ = LinearInterpolation((θ_default, ϕ_default), Gθ),
        ϕ = LinearInterpolation((θ_default, ϕ_default), Gϕ)
    )
end

# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(point::Vector{anten_point}, point_I, θₜ, ϕₜ, k = k, θ = θ_default, ϕ = ϕ_default)

    # # apply_rotation to the pattern
    # @time map(
    #     (point::anten_point) ->
    #         point.pattern = rotate_pattern(point.pattern, point.local_coord),
    #     point
    # )
    @time Threads.@threads for p in point 
        p.pattern = rotate_pattern(p.pattern, p.local_coord)
    end


    # calculate result
    result_θ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1), size(point, 1))
    result_ϕ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1), size(point, 1))
    Threads.@threads for (ind_point, point_i) = collect(enumerate(point))
        for (i, θ) = enumerate(θ), (j, ϕ) = enumerate(ϕ)
            result_θ[i, j, ind_point] = Iₛ(point_I[ind_point].p, θₜ, ϕₜ, k)Pi(point_i.p, point_i.pattern.θ, θ, ϕ, θₜ, ϕₜ, k)
            result_ϕ[i, j, ind_point] = Iₛ(point_I[ind_point].p, θₜ, ϕₜ, k)Pi(point_i.p, point_i.pattern.ϕ, θ, ϕ, θₜ, ϕₜ, k)
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
function anten_read(filepath, type = "hfss"; unit = "abs", factor = 1, order = [2, 1, 4, 3])
    #pick data from the dataframe
    unit_corrector = unit -> begin
        if unit == "db"
            x -> 10^(x / 10) |>sqrt
        else
            x -> x * factor
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
# watt/unit solid angle 2-12a
radiation_intensity = pattern -> [
    1/2(120pi)*(abs(pattern.θ(θ, ϕ))^2 + abs(pattern.ϕ(θ, ϕ))^2)
    for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
# watt 2-13
radiated_power = pattern -> P = [
    sin(θ)*U *2pi *pi  * 1/size(ϕ_default,1)size(θ_default, 1)
    for (U, θ) in zip(radiation_intensity(pattern),θ_grid) ] |> sum
# 2-16
directivity = pattern ->  4pi*radiation_intensity(pattern)/radiated_power(pattern) 
gain = (pattern, inciden_power) -> 4pi*radiation_intensity(pattern)/inciden_power
end
