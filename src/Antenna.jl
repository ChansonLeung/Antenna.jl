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
using Peaks

using LoopVectorization
using Tullio

include("type.jl")
include("utils.jl")
include("plotting.jl")
include("read_re_file.jl")

export
    cal_pattern,
    anten_read,
    rotate_vec_in_sph,
    rotate_vec_in_sph!,
    directivity,
    radiation_intensity,
    radiated_power,
    SLL_maxgain,
    gain,
    Iₛ,
    vec_θ,
    vec_θ_static,
    vec_ϕ,
    vec_ϕ_static,
    #degbug
    rotate_pattern_tullion

# array function
# array factor
AF(point, θ, ϕ, k) = exp(1im * k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
# current source
Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k)^-1
# Iₛ(point, θₜ, ϕₜ, k) = 1
# apply rotation
# XXX untest
function rotate_vec_in_sph(θ, ϕ, M)
    point = sph2cart(θ, ϕ, 1)
    vec = M * point
    θ1, ϕ1, r = cart2sph(vec[1], vec[2], vec[3])
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

function rotate_vec_in_sph!(point, M)
    sph2cart!(point, point[1], point[2], 1)
    mul!(point, M, @SArray [point[1], point[2], point[3]])
    cart2sph!(point, point[1], point[2], point[3])
end

vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
vec_θ_static(θ, ϕ) = @SArray [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]

vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0.0]
vec_ϕ_static(θ, ϕ) = @SArray [-sin(ϕ), cos(ϕ), 0.0]

function cart2sph_static(x, y, z)
    # XXX r are not calcuated, may cause problem
    r = sqrt(x^2 + y^2 + z^2)
    θ = acos(z / r)
    ϕ = atan(y, x)
    @SVector[θ, ϕ, r]
end

function sph2cart_static(θ, ϕ, r)
    @SArray [r * sin(θ)cos(ϕ),
        r * sin(θ)sin(ϕ),
        r * cos(θ)]
end

function rotate_pattern_tullion(res, pattern, coord, interp_exported=true)
    res = reinterpret(reshape, SVector{2,ComplexF64}, res)
    mat = SMatrix{3,3}(coord)
    patternθ = pattern.θ
    patternϕ = pattern.ϕ
    θ = θ_default
    ϕ = ϕ_default
    # close threads is needed, otherwise lead to segmentation fault(may due to this unformally supported sytle for @tullio)
    @tullio threads = false res[i, j] = begin
        # rotate grid
        sinθ, sinϕ, cosθ, cosϕ = sin(θ[i]), sin(ϕ[j]), cos(θ[i]), cos(ϕ[j])
        x, y, z = SVector{3}(mat' * @SVector[sinθ * cosϕ, sinθ * sinϕ, cosθ])
        θ_rot, ϕ_rot, r_rot = cart2sph_static(x, y, z)
        # get pattern vector value in rotate grid
        sinθ_rot, sinϕ_rot, cosθ_rot, cosϕ_rot = sin(θ_rot), sin(ϕ_rot), cos(θ_rot), cos(ϕ_rot)
        val_vec_θ_rot = SVector{3}(mat * patternθ(θ_rot, ϕ_rot) * @SVector[cosθ_rot * cosϕ_rot, cosθ_rot * sinϕ_rot, -sinθ_rot])
        val_vec_ϕ_rot = SVector{3}(mat * patternϕ(θ_rot, ϕ_rot) * @SVector[-sin(ϕ_rot), cos(ϕ_rot), 0.0])
        # project vector to global grid
        vec_θ_global = @SVector[cosθ * cosϕ, cosθ * sinϕ, -sinθ]
        vec_ϕ_global = @SVector[-sin(ϕ[j]), cos(ϕ[j]), 0.0]
        @SVector[dot(val_vec_θ_rot, vec_θ_global) + dot(val_vec_ϕ_rot, vec_θ_global),
            dot(val_vec_ϕ_rot, vec_ϕ_global) + dot(val_vec_θ_rot, vec_ϕ_global)]
    end
    if interp_exported
        res = reinterpret(reshape, ComplexF64, res)
        anten_pattern(
            θ=linear_interpolation((θ, ϕ), @view res[1, :, :]),
            ϕ=linear_interpolation((θ, ϕ), @view res[2, :, :])
        )
    end
end


# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(point::Vector{anten_point}, θₜ::Float64, ϕₜ::Float64; k::Float64=k, θ=θ_default, ϕ=ϕ_default, spin=false)
    spin && @floop for p in point
        # p.pattern_grid = rotate_pattern_tullion(p.pattern_grid, p.pattern, p.local_coord)
        rotate_pattern_tullion(p.pattern_grid, p.pattern, p.local_coord)
    end

    positions = getfield.(point, :p)
    coeffi = getfield.(point, :coeffi)
    pattern = getfield.(point, :pattern_grid)

    @tullio I[p] := Iₛ(positions[p], θₜ, ϕₜ, k)
    @tullio result[i, j] := coeffi[p] * I[p] * AF(positions[p], θ[i], ϕ[j], k) * SVector{2}((@view pattern[p][:, i, j]))
    result = reinterpret(reshape, ComplexF64, result)
    anten_pattern(
        θ=linear_interpolation((θ, ϕ), @view result[1, :, :]),
        ϕ=linear_interpolation((θ, ϕ), @view result[2, :, :])
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


function SLL_maxgain(pattern::anten_pattern; θ=0.0, ϕ=0.0)
    d = directivity(pattern)
    curve = [d((mod_angle_deg(θ, ϕ) .|> deg2rad)...) for θ in -100:100, ϕ in ϕ]
    val = maxima(curve) |> sort!
    val[end] / val[end-1], val[end]
end

end
