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
using CUDA

include("type.jl")
include("utils.jl")
include("plotting.jl")
include("read_re_file.jl")
include("makie_antenna.jl")
include("typical_shape.jl")

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
    maximum_directivity,
    linear_interp_pattern,
    linear_interp_pattern_2com,
    #degbug
    rotate_pattern_tullion,
    cal_pattern_cuda,
    radiated_power_beta,
    directivity_beta,
    radiation_intensity_beta

# array function
# array factor
# AF(point, θ, ϕ, k) = exp(1im * k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
AF(point, θ, ϕ, k) = cis(k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
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

function rotate_pattern_tullion(res, pattern, coord; interp=true)
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
        if θ_rot ≈ 0
        # θ = 45, ϕ = 90
        # θ_rot = 0, ϕ = ? 
            ϕ_rot =ϕ[j] 
        end

        # get pattern vector value in rotate grid
        sinθ_rot, sinϕ_rot, cosθ_rot, cosϕ_rot = sin(θ_rot), sin(ϕ_rot), cos(θ_rot), cos(ϕ_rot)
        val_vec_θ_rot = SVector{3}(mat * patternθ(θ_rot, ϕ_rot) * @SVector[cosθ_rot * cosϕ_rot, cosθ_rot * sinϕ_rot, -sinθ_rot])
        val_vec_ϕ_rot = SVector{3}(mat * patternϕ(θ_rot, ϕ_rot) * @SVector[-sinϕ_rot, cosϕ_rot, 0.0])
        # project vector to global grid
        vec_θ_global = @SVector[cosθ * cosϕ, cosθ * sinϕ, -sinθ]
        vec_ϕ_global = @SVector[-sin(ϕ[j]), cos(ϕ[j]), 0.0]

        @SVector[dot(conj.(val_vec_θ_rot), vec_θ_global) + dot(conj.(val_vec_ϕ_rot), vec_θ_global),
            dot(conj.(val_vec_ϕ_rot), vec_ϕ_global) + dot(conj.(val_vec_θ_rot), vec_ϕ_global)]
    end
    if interp
        res = reinterpret(reshape, ComplexF64, res)
        anten_pattern(
            θ=interpolate!((θ, ϕ), (@view res[1, :, :]), Gridded(Linear())),
            ϕ=interpolate!((θ, ϕ), (@view res[2, :, :]), Gridded(Linear()))
        )
    end
end


# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(
        point::Vector{anten_point},
        θₜ,
        ϕₜ;
        k::Float64=k,
        θ=θ_default,
        ϕ=ϕ_default,
        spin=false,
        I=nothing,
        withAF=true,
        optimal_directivity=nothing,
        result = nothing,
        args...
    )
    θₜ = Float64(θₜ)
    ϕₜ = Float64(ϕₜ)
    D_direct = ones(length(point))

    spin &&  for (idx,p) in enumerate(point)
        pattern_i = if p.pattern_grid_cache == false
            p.pattern_grid_cache = true
            pattern_i = rotate_pattern_tullion(p.pattern_grid, p.pattern, p.local_coord)
        else
            linear_interp_pattern_2com(p.pattern_grid)
        end
        (optimal_directivity !== nothing) &&  (
            D_direct[idx] = directivity_beta(
                pattern_i, 
                θₜ, 
                ϕₜ, 
                optimal_directivity, power=radiated_power_beta(p.pattern)
            )
        )
    end


    positions = getfield.(point, :p)
    coeffi = getfield.(point, :coeffi) .* sqrt.(D_direct)
    pattern = getfield.(point, :pattern_grid)
    pattern_with_AF = getfield.(point, :pattern_grid_withAF)

    I===nothing && @tullio I[p] := Iₛ(positions[p], θₜ, ϕₜ, k)

    if withAF
      @floop  for idx_point in  eachindex(pattern_with_AF)
            if point[idx_point].pattern_grid_withAF_cache==false
                pattern_with_AF_i = pattern_with_AF[idx_point]
                pattern_i = pattern[idx_point]
                for idx_grid = CartesianIndices(θ_grid)
                    tmp_AF = AF(positions[idx_point], θ_grid[idx_grid], ϕ_grid[idx_grid], k)
                    pattern_with_AF_i[1, idx_grid[1],idx_grid[2]] = tmp_AF * pattern_i[1, idx_grid[1],idx_grid[2]]
                    pattern_with_AF_i[2, idx_grid[1],idx_grid[2]] = tmp_AF * pattern_i[2, idx_grid[1],idx_grid[2]]
                end
                point[idx_point].pattern_grid_withAF_cache = true
            end
        end
    end


    if isnothing(result)
        result = zeros(SVector{2, ComplexF64},size(θ_grid)...)
    end
    # @tullio result[i, j] = coeffi[p] * I[p] * AF(positions[p], θ[i], ϕ[j], k) * SVector{2}((@view pattern[p][:, i, j]))
    @tullio result[i, j] = coeffi[p] * I[p] * SVector{2}((@view pattern_with_AF[p][:, i, j]))
    result = reinterpret(reshape, ComplexF64, result)

    return anten_pattern(
        θ=interpolate!((θ, ϕ), (@view result[1, :, :]), Gridded(Linear())),
        ϕ=interpolate!((θ, ϕ), (@view result[2, :, :]), Gridded(Linear()))
    ), I.*coeffi
end

linear_interp_pattern(pattern_grid) = interpolate!((θ_default, ϕ_default), pattern_grid, Gridded(Linear()))
linear_interp_pattern_2com(pattern_grid) = begin
    result = reinterpret(reshape, ComplexF64, pattern_grid)
    return anten_pattern(
        θ=interpolate!((θ_default, ϕ_default), (@view result[1, :, :]), Gridded(Linear())),
        ϕ=interpolate!((θ_default, ϕ_default), (@view result[2, :, :]), Gridded(Linear()))
    )
end


function cal_pattern_cuda_kernel(
    positions, 
    coeffi, 
    mat_rot, 
    pattern_in_re, 
    pattern_in_im, 
    pattern_out, 
    k,
    θₜ=0, 
    ϕₜ=0, 
    θ_default=cu(LinRange(0, pi, 181)), 
    ϕ_default=cu(LinRange(-pi, pi, 361)),
)
    idx_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride_i = gridDim().x * blockDim().x

    for i in idx_i:stride_i:length(θ_default)*length(ϕ_default)
        for j in 1:length(coeffi)
            mat_rot_i = SMatrix{3,3}(mat_rot[j])
            idx_θ = (i - 1) % length(θ_default) + 1
            idx_ϕ = (i - 1) ÷ length(θ_default) + 1

            # rotate grid
            θ, ϕ = θ_default[idx_θ], ϕ_default[idx_ϕ]
            sinθ, sinϕ, cosθ, cosϕ = sin(θ), sin(ϕ), cos(θ), cos(ϕ)
            x, y, z = mat_rot_i' * @SVector[sinθ * cosϕ, sinθ * sinϕ, cosθ]
            θ_rot, ϕ_rot, r_rot = cart2sph_static(x, y, z)

            # get pattern vector value in rotate grid
            sinθ_rot, sinϕ_rot, cosθ_rot, cosϕ_rot = sin(θ_rot), sin(ϕ_rot), cos(θ_rot), cos(ϕ_rot)
                # FIXME need interpolation
            tex_idx_θ = θ_rot * (length(θ_default)-1)/pi+1
            tex_idx_ϕ = (ϕ_rot +pi) * (length(ϕ_default)-1)/2pi+1
            Gθ_rot_local = pattern_in_re[1,tex_idx_θ, tex_idx_ϕ] + 1im*pattern_in_im[1,tex_idx_θ, tex_idx_ϕ]
            Gϕ_rot_local = pattern_in_re[2,tex_idx_θ, tex_idx_ϕ] + 1im*pattern_in_im[2,tex_idx_θ, tex_idx_ϕ]

            val_vec_θ_rot =  Gθ_rot_local* mat_rot_i * @SVector[cosθ_rot * cosϕ_rot, cosθ_rot * sinϕ_rot, -sinθ_rot]
            val_vec_ϕ_rot =  Gϕ_rot_local* mat_rot_i * @SVector[-sinϕ_rot, cosϕ_rot, 0.0]

            # project vector to global grid
            vec_θ_global = @SVector[cosθ * cosϕ, cosθ * sinϕ, -sinθ]
            vec_ϕ_global = @SVector[-sinϕ, cosϕ, 0.0]

            G_θ = dot(val_vec_θ_rot, vec_θ_global) + dot(val_vec_ϕ_rot, vec_θ_global)
            G_ϕ = dot(val_vec_ϕ_rot, vec_ϕ_global) + dot(val_vec_θ_rot, vec_ϕ_global)

            sinθₜ, sinϕₜ, cosθₜ, cosϕₜ = sin(θₜ), sin(ϕₜ), cos(θₜ), cos(ϕₜ)
            p = @SVector[positions[1, j], positions[2, j], positions[3, j]]
            # # ---
            # tmp = coeffi[j] * Iₛ(p, θₜ, ϕₜ, k) * AF(p, θ, ϕ, k)
            tmp = coeffi[j] * exp(1.0im * k * (p[1] * sinθ*cosϕ + p[2] * sinθ*sinϕ + p[3] * cosθ))
            # exp(1.0im * k * (p[1] * sinθₜ*cosϕₜ + p[2] * sinθₜ*sinϕₜ + p[3] * cosθₜ))

            pattern_out[1, idx_θ, idx_ϕ] += tmp * G_θ
            pattern_out[2, idx_θ, idx_ϕ] += tmp * G_ϕ
        end
    end
    return nothing
end


function cal_pattern_cuda(point::Vector{anten_point}, k::Float64=k, θ=θ_default, ϕ=ϕ_default, spin=false)
    positions = reduce(hcat,getfield.(point,:p))
    coeffi = reduce(vcat,getfield.(point,:coeffi))
    mat_rot = [SMatrix{3,3}(i.local_coord) for i in point]
    pattern_in_re = CuTextureArray(cu(real.(point[1].pattern_grid))) |> x -> CuTexture(x; interpolation=CUDA.LinearInterpolation())
    pattern_in_im = CuTextureArray(cu(imag.(point[1].pattern_grid))) |> x -> CuTexture(x; interpolation=CUDA.LinearInterpolation())
    pattern_out_cu = zeros(ComplexF32, 2, size(θ_grid)...) |> cu

    @cuda threads = 128 blocks=ceil(Int, length(θ_default)*length(ϕ_default)/128) cal_pattern_cuda_kernel(
        cu(positions),
        cu(coeffi),
        cu(mat_rot),
        pattern_in_re,
        pattern_in_im,
        pattern_out_cu,
        k
);
    result = Array(pattern_out_cu)

    anten_pattern(
        θ=interpolate!((θ, ϕ), (@view result[1, :, :]), Gridded(Linear())),
        ϕ=interpolate!((θ, ϕ), (@view result[2, :, :]), Gridded(Linear()))
    )
end

# function below may has performance problem that need to be check because of the usage of lambda inside the comprehension 
# unit of radiation intensity is w/sr (watt/unit solid angle) refer to Antenna Theory 2-12a
# FIXME too slow!!
function radiation_intensity(pattern::anten_pattern, component=:all)
    if component == :all
        (θ, ϕ) -> 1 / 2(120pi) * (abs(pattern.θ(θ, ϕ))^2 + abs(pattern.ϕ(θ, ϕ))^2)
    elseif component==:θ
        (θ, ϕ) -> 1 / 2(120pi) * (abs(pattern.θ(θ, ϕ))^2)
    elseif component==:ϕ
        (θ, ϕ) -> 1 / 2(120pi) * (abs(pattern.ϕ(θ, ϕ))^2)
    end
end
function radiation_intensity_beta(pattern::anten_pattern, θ, ϕ, component=:all)
    patternθ = pattern.θ
    patternϕ = pattern.ϕ
    if length(θ)>1 && length(ϕ)>1
        if component == :all
            @tullio res[i,j] := 1 / 2(120pi) * (abs(patternθ(θ[i,j], ϕ[i,j]))^2 + abs(patternϕ(θ[i,j], ϕ[i,j]))^2)
        elseif component==:θ
            @tullio res[i,j] := 1 / 2(120pi) * abs(patternθ(θ[i,j], ϕ[i,j]))^2 
        elseif component==:ϕ
            @tullio res[i,j] := 1 / 2(120pi) * abs(patternϕ(θ[i,j], ϕ[i,j]))^2
        elseif component==:RHCP
            @tullio res[i,j] := begin
                θi = θ[i,j]
                ϕi = ϕ[i,j]
                vec_θ = vec_θ_static(θi,ϕi)
                vec_ϕ = vec_ϕ_static(θi,ϕi)
                E_vec = patternθ(θi,ϕi)*vec_θ +  patternϕ(θi,ϕi)*vec_ϕ
                RHCP_base = (vec_θ_static(θi,ϕi) - 1im*vec_ϕ_static(θi,ϕi))./sqrt(2)
                1/2(120pi) * (dot(E_vec,RHCP_base))
            end
        end
    else
        if component == :all
            (@. 1 / 2(120pi) * (abs(patternθ(θ, ϕ))^2 + abs(patternϕ(θ, ϕ))^2))
        elseif component==:θ
            @. 1 / 2(120pi) * abs(patternθ(θ, ϕ))^2 
        elseif component==:ϕ
            @. 1 / 2(120pi) * abs(patternϕ(θ, ϕ))^2
        end
    end
end
# unit of radiated power is w (watt) refer to Antenna Theory 2-13
@memoize function radiated_power(pattern::anten_pattern, component=:all)
    result = sum([sin(θ) * U * 2pi * pi * 1 / size(ϕ_default, 1)size(θ_default, 1)
                  for (U, θ) in zip(radiation_intensity(pattern, component).(θ_grid, ϕ_grid), θ_grid)])
end
#FIXME the radiated_power_beta did not take component into account
function radiated_power_beta(pattern::anten_pattern, component=:all)
    if pattern.radiated_power === nothing
        pattern.radiated_power = sum([sin(θ) * U * 2pi * pi * 1 / size(ϕ_default, 1)size(θ_default, 1)
                  for (U, θ) in zip(radiation_intensity(pattern, component).(θ_grid, ϕ_grid), θ_grid)])
    else
        pattern.radiated_power
    end
end

# no unit refer to Antenna Theory 2-16
@memoize function directivity(pattern::anten_pattern, component=:all)
    (θ, ϕ) -> 4pi * radiation_intensity(pattern, component)(θ, ϕ) / radiated_power(pattern)
end
function directivity_beta(pattern::anten_pattern, θ,ϕ, component=:all; power=nothing)
    if power === nothing
       4pi * radiation_intensity_beta(pattern, θ, ϕ, component) / radiated_power_beta(pattern)
    else
       4pi * radiation_intensity_beta(pattern, θ, ϕ, component) / power
    end
end
function maximum_directivity(pattern::anten_pattern)
   dir = directivity(pattern).(θ_grid, ϕ_grid)
   max = maximum(dir)
   max, (θ_grid[dir.==max][1] , ϕ_grid[dir.==max][1])
end


# no unit refer to Antenna Theory 2-16
gain = (pattern, inciden_power) -> function (θ, ϕ, component=:all)
    4pi * radiation_intensity(pattern, component)(θ, ϕ) / inciden_power
end


function SLL_maxgain(pattern::anten_pattern; θ=0.0, ϕ=0.0)
    d = directivity(pattern)
    curve = [d((mod_angle_deg(θ, ϕ) .|> deg2rad)...) for θ in -100:100, ϕ in ϕ]
    val = maxima(curve) |> sort!
    val[end] / val[end-1], val[end]
end

end
