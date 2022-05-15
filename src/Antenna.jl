module Antenna

using LinearAlgebra
using Interpolations
using Chain
# pattern file
using CSV
using DataFrames
using Match
using Reexport
using Memoize
using FLoops

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
AF(point, θ, ϕ, k) = exp(1im * k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
# current source
Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k)^-1
# Iₛ(point, θₜ, ϕₜ, k) = 1
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
rotate_vec_in_cart = (vec, M) -> M * vec
vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0]


rotate_pattern = (pattern::anten_pattern, coord::Matrix{Float64}) -> begin

    θ_grid::Matrix{Float64} = [θ for θ in θ_default, ϕ in ϕ_default]
    ϕ_grid::Matrix{Float64} = [ϕ for θ in θ_default, ϕ in ϕ_default]

    Gθ_grid::Matrix{Vector{ComplexF64}} = [
        begin
            θ, ϕ = rotate_vec_in_sph(θ, ϕ, coord')
            (coord * vec_θ(θ, ϕ)) * pattern.θ(θ, ϕ)
        end
        for (θ, ϕ) = zip(θ_grid, ϕ_grid)
    ]

    Gϕ_grid::Matrix{Vector{ComplexF64}} = [
        begin
            θ, ϕ = rotate_vec_in_sph(θ, ϕ, coord')
            (coord * vec_ϕ(θ, ϕ)) * pattern.ϕ(θ, ϕ)
        end
        for (θ, ϕ) = zip(θ_grid, ϕ_grid)
    ]

    vec_θ₁_map_grid::Matrix{Vector{Float64}} = [vec_θ(θ, ϕ) for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    vec_ϕ₁_map_grid::Matrix{Vector{Float64}} = [vec_ϕ(θ, ϕ) for (θ, ϕ) = zip(θ_grid, ϕ_grid)]

    Gθ = (dot.(Gθ_grid, vec_θ₁_map_grid)) .+
         (dot.(Gϕ_grid, vec_θ₁_map_grid))
    Gϕ = (dot.(Gϕ_grid, vec_ϕ₁_map_grid)) .+
         (dot.(Gθ_grid, vec_ϕ₁_map_grid))

    result = anten_pattern(
        θ=LinearInterpolation((θ_default, ϕ_default), Gθ),
        ϕ=LinearInterpolation((θ_default, ϕ_default), Gϕ)
    )
end

# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(point::Vector{anten_point},  θₜ, ϕₜ, k=k, θ=θ_default, ϕ=ϕ_default)

    # # apply_rotation to the pattern
    # @time map(
    #     (point::anten_point) ->
    #         point.pattern = rotate_pattern(point.pattern, point.local_coord),
    #     point
    # )
    # @time "rotate" Threads.@threads for p in point 
    #     p.pattern = rotate_pattern(p.pattern, p.local_coord)
    # end

    # calculate result
    result_θ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1))
    result_ϕ = zeros(ComplexF64, size(θ, 1), size(ϕ, 1))
    @floop for p_i in point
        pattern_θ = p_i.pattern.θ
        pattern_ϕ = p_i.pattern.ϕ
        tmp = p_i.coeffi * Iₛ(p_i.p, θₜ, ϕₜ, k) .* [AF(p_i.p, θi, ϕi, k)  for θi in θ_default, ϕi in ϕ_default]
      @reduce result_θ .+=  tmp.* [pattern_θ(θi, ϕi)  for θi in θ_default, ϕi in ϕ_default]
      @reduce result_ϕ .+=  tmp.* [pattern_ϕ(θi, ϕi)  for θi in θ_default, ϕi in ϕ_default]
    end

    anten_pattern(
        θ=LinearInterpolation((θ_default, ϕ_default), result_θ .|>abs),
        ϕ=LinearInterpolation((θ_default, ϕ_default), result_ϕ .|>abs),
    )
end

# order: the column number for θ,ϕ,Sθ,Sϕ, for HFSS is [2, 1, 4, 3]
function anten_read(filepath, type="hfss"; unit="abs", factor=1, order=[2, 1, 4, 3])
    #pick data from the dataframe
    unit_corrector = unit -> begin
        if unit == "db"
            x -> 10^(x / 10) |> sqrt
        else
            x -> x * factor
        end
    end
    convert_from_hfss = x -> begin
        (
            θ=deg2rad.(x[:, order[1]]),
            ϕ=deg2rad.(x[:, order[2]]),
            Gθ=x[:, order[3]] .|> unit_corrector(unit),
            Gϕ=x[:, order[4]] .|> unit_corrector(unit)
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
        θ=translate2Interpolation(θ, ϕ, raw_anten.Gθ),
        ϕ=translate2Interpolation(θ, ϕ, raw_anten.Gϕ)
    )
end

# function below may has performance problem that need to be check because of the usage of lambda inside the comprehension 
# unit of radiation intensity is w/sr (watt/unit solid angle) refer to Antenna Theory 2-12a
function radiation_intensity(pattern::anten_pattern)
    (θ, ϕ) -> 1 / 2(120pi) * (abs(pattern.θ(θ, ϕ))^2 + abs(pattern.ϕ(θ, ϕ))^2)
end
# unit of radiated power is w (watt) refer to Antenna Theory 2-13
@memoize function radiated_power(pattern::anten_pattern)
    result = sum([sin(θ) * U * 2pi * pi * 1 / size(ϕ_default, 1)size(θ_default, 1)
                  for (U, θ) in zip(radiation_intensity(pattern).(θ_grid, ϕ_grid), θ_grid)])
end

# no unit refer to Antenna Theory 2-16
@memoize function directivity(pattern::anten_pattern)
    (θ, ϕ) -> 4pi * radiation_intensity(pattern)(θ, ϕ) / radiated_power(pattern)
end
# no unit refer to Antenna Theory 2-16
gain = (pattern, inciden_power) -> function (θ, ϕ)
    4pi * radiation_intensity(pattern)(θ, ϕ) / inciden_power
end
end
