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

include("type.jl")
include("utils.jl")
include("plotting.jl")
include("read_re_file.jl")

export
    cal_pattern,
    anten_read,
    rotate_vec_in_sph,
    rotate_vec_in_sph!,
    rotate_pattern,
    rotate_pattern!,
    directivity,
    radiation_intensity,
    radiated_power,
    SLL_maxgain,
    gain,
    Iₛ,
    vec_θ,
    vec_θ_static,
    vec_ϕ,
    vec_ϕ_static

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
    cart2sph!(point, point[1],point[2], point[3])
end

vec_θ(θ, ϕ) = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
vec_θ_static(θ, ϕ) = @SArray [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
# vec_θ(θ, ϕ) = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), -sin(θ)]
# vec_θ_static(θ, ϕ) = @SArray [sin(θ)cos(ϕ), sin(θ)sin(ϕ), -sin(θ)]

vec_ϕ(θ, ϕ) = [-sin(ϕ), cos(ϕ), 0.0]
vec_ϕ_static(θ, ϕ) = @SArray [-sin(ϕ), cos(ϕ), 0.0]


function rotate_pattern(pattern::anten_pattern, coord::Matrix{Float64})

    # θ_grid::Matrix{Float64} = [θ for θ in θ_default, ϕ in ϕ_default]
    # ϕ_grid::Matrix{Float64} = [ϕ for θ in θ_default, ϕ in ϕ_default]
    
    inverse_rotate_grid = rotate_vec_in_sph.(θ_grid, ϕ_grid, [coord'])
    inv_rot_θ = map(x->x[1], inverse_rotate_grid)
    inv_rot_ϕ = map(x->x[2], inverse_rotate_grid)

    Gθ_grid::Matrix{Vector{ComplexF64}} = [
        begin
            (coord*vec_θ(θ, ϕ)) * pattern.θ(θ, ϕ)
        end
        for (θ, ϕ) = zip(inv_rot_θ, inv_rot_ϕ)
    ]

    Gϕ_grid::Matrix{Vector{ComplexF64}} = [
        begin
            (coord*vec_ϕ(θ, ϕ)) * pattern.ϕ(θ, ϕ)
        end
        for (θ, ϕ) = zip(inv_rot_θ, inv_rot_ϕ)
    ]

    vec_θ₁_map_grid::Matrix{Vector{Float64}} = [vec_θ(θ, ϕ) for (θ, ϕ) = zip(θ_grid, ϕ_grid)]
    vec_ϕ₁_map_grid::Matrix{Vector{Float64}} = [vec_ϕ(θ, ϕ) for (θ, ϕ) = zip(θ_grid, ϕ_grid)]

    Gθ = (dot.(Gθ_grid, vec_θ₁_map_grid)) .+
         (dot.(Gϕ_grid, vec_θ₁_map_grid))
    Gϕ = (dot.(Gϕ_grid, vec_ϕ₁_map_grid)) .+
         (dot.(Gθ_grid, vec_ϕ₁_map_grid))

    result = anten_pattern(
        θ=linear_interpolation((θ_default, ϕ_default), Gθ, extrapolation_bc = Flat()),
        ϕ=linear_interpolation((θ_default, ϕ_default), Gϕ, extrapolation_bc = Flat())
    )
end
function rotate_pattern!(pattern::anten_pattern, coord::Matrix{Float64})
    # convert 3xMxN matrix to vector{eltpe, (M*N)}
    point_mode(grid) = reshape(grid, 3,:)|>eachcol
    vec2grid(vec) = reshape(vec, size(ϕ_grid)...)

    inv_rot_grid = hcat(vec(θ_grid), vec(ϕ_grid), ones(length(ϕ_grid)))
    for row in eachrow(inv_rot_grid)
        rotate_vec_in_sph!(row, coord')
    end

    inv_rot_θ = @view inv_rot_grid[:,1]
    inv_rot_ϕ = @view inv_rot_grid[:,2]
    inv_rot_θ = reshape(inv_rot_θ, size(θ_grid))
    inv_rot_ϕ = reshape(inv_rot_ϕ, size(ϕ_grid))

    Gainθ_grid = pattern.θ.(inv_rot_θ,inv_rot_ϕ)
    vec_Gainθ_grid = zeros(3,size(inv_rot_θ)...)
    for (point,θ, ϕ) in zip(point_mode(vec_Gainθ_grid), vec(inv_rot_θ), vec(inv_rot_ϕ))
        point .= vec_θ_static(θ,ϕ)
        mul!(point, coord, @SArray[point[1], point[2], point[3]])
    end
    Gainθ_grid = vec_Gainθ_grid .* reshape(Gainθ_grid,1,size(Gainθ_grid)... )


    Gainϕ_grid = pattern.ϕ.(inv_rot_θ,inv_rot_ϕ)
    #preallocate memory
    vec_Gainϕ_grid = zeros(3,size(inv_rot_ϕ)...)
    for (point,θ, ϕ) in zip(point_mode(vec_Gainϕ_grid), vec(inv_rot_θ), vec(inv_rot_ϕ))
        point .= vec_ϕ_static(θ,ϕ)
        mul!(point, coord, @SArray[point[1], point[2], point[3]])
    end
    Gainϕ_grid = vec_Gainϕ_grid  .*reshape(Gainϕ_grid,1,size(Gainϕ_grid)...)

    #preallocate memory
    vec_θ₁_map_grid = zeros(3,size(θ_grid)...)
    vec_ϕ₁_map_grid = zeros(3,size(θ_grid)...)
    for (point_θ, point_ϕ, θ, ϕ) in zip(point_mode(vec_θ₁_map_grid), point_mode(vec_ϕ₁_map_grid),θ_grid, ϕ_grid)
        point_θ .= vec_θ_static(θ,ϕ)
        point_ϕ .= vec_ϕ_static(θ,ϕ)
    end

    prealloc_to_grid = vec2grid ∘ collect ∘ point_mode
    Gθ = (dot.(Gainθ_grid |>prealloc_to_grid, vec_θ₁_map_grid|>prealloc_to_grid)) .+
         (dot.(Gainϕ_grid |>prealloc_to_grid, vec_θ₁_map_grid|>prealloc_to_grid))
    Gϕ = (dot.(Gainϕ_grid |>prealloc_to_grid, vec_ϕ₁_map_grid|>prealloc_to_grid)) .+
         (dot.(Gainθ_grid |>prealloc_to_grid, vec_ϕ₁_map_grid|>prealloc_to_grid))

    result = anten_pattern(
        θ=linear_interpolation((θ_default, ϕ_default), Gθ, extrapolation_bc = Flat()),
        ϕ=linear_interpolation((θ_default, ϕ_default), Gϕ, extrapolation_bc = Flat())
    )
end

# calculate the global pattern
# para∑ can be write like this
# para∑(Pi(point:p, point:pattern.θ , θ, ϕ, θₜ, ϕₜ, k), (point, θ,ϕ, θₜ, ϕₜ,k))
function cal_pattern(point::Vector{anten_point},  θₜ::Float64, ϕₜ::Float64; k::Float64=k , θ=θ_default, ϕ=ϕ_default, spin=false)

    # apply_rotation to the pattern
    # point = []
    # @floop for p in point_in
    #     push!(point,  anten_point(p = p.p, pattern =rotate_pattern!(p.pattern, p.local_coord)))
    # end

    spin && @floop for p in point
        p.pattern = rotate_pattern!(p.pattern, p.local_coord)
    end



    # calculate result
    result_θ = zeros(ComplexF64, length(θ), length(ϕ))
    result_ϕ = zeros(ComplexF64, length(θ), length(ϕ))
    tmp =  zeros(ComplexF64, length(θ), length(ϕ))
    
    for p_i in point
        pattern_θ = p_i.pattern.θ
        pattern_ϕ = p_i.pattern.ϕ
        co = p_i.coeffi * Iₛ(p_i.p, θₜ, ϕₜ, k)
       @inbounds for index_c in  CartesianIndices(tmp)
            θi = θ[index_c[1]]
            ϕi = ϕ[index_c[2]]
            tmp[index_c] = co * AF(p_i.p, θi, ϕi, k)
        end
        result_θ .+= tmp.* pattern_θ.(θ_grid,ϕ_grid)
        result_ϕ .+= tmp.* pattern_ϕ.(θ_grid,ϕ_grid)
    end
    
    anten_pattern(
        θ=linear_interpolation((θ, ϕ), result_θ .|>abs, extrapolation_bc = Flat()) ,
        ϕ=linear_interpolation((θ, ϕ), result_ϕ .|>abs, extrapolation_bc = Flat()) 
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


function SLL_maxgain(pattern::anten_pattern;θ=0.0, ϕ=0.0)
    d = directivity(pattern)
    curve = [d((mod_angle_deg(θ,ϕ) .|>deg2rad)...) for θ in -100:100, ϕ in ϕ]
    val = maxima(curve) |> sort! 
    val[end]/val[end-1], val[end]
end

end
