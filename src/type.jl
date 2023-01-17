# FIXME
using LinearAlgebra
include("synthesis.jl")

# XXX only for uniform gird
# 1. read the pattern csv file and plot it
# elem_pattern = anten_read("pattern.csv")
# plot(pattern)
#
# 2. calculate the interpolated
# pattern = [elem_pattern.pattern(θ,ϕ) for θ=θ_default, ϕ=ϕ_default]
Base.@kwdef mutable struct anten_pattern
    θ::Any
    ϕ::Any
end

# # XXX the linear Interpolations may access the extrapolate point because of the accuracy of float numer
# @recipe ploting(p::anten_pattern; ϕ::Float64 = pi/2-0.00001) = begin
#     squre_sum_sqrt = (x, y) -> (sqrt(abs(x)^2 + abs(y)^2) |> x -> 20log10(x))
#     # squre_sum_sqrt = (x, y) -> 20log10(x)
#     (
#         [θ_default -θ_default] .|>rad2deg,
#         [
#             [squre_sum_sqrt(p.θ(θ, ϕ), p.ϕ(θ, ϕ)) for θ = θ_default, ϕ = ϕ];;
#             [squre_sum_sqrt(p.θ(θ, ϕ), p.ϕ(θ, ϕ)) for θ = θ_default, ϕ = mod2pi(ϕ - pi + pi) -pi]
#         ]
#     )
# end
# pattern_identity = anten_pattern(θ = (θ, ϕ) -> abs(cos(θ)^2), ϕ = (θ, ϕ) -> 0)
η = 120pi
pattern_identity = anten_pattern(θ=(θ, ϕ) -> 1, ϕ=(θ, ϕ) -> 0)
# antenna theory 4-84
pattern_dipole = anten_pattern(
    θ=(θ, ϕ) -> 1im * η * exp(-1im * k) / (2pi) * (cos(pi / 2 * cos(θ)) / (sin(θ) + 1e-6)),
    ϕ=(θ, ϕ) -> 0
)

Base.@kwdef mutable struct anten_point 
    #point
    p::Vector{Float64}
    pattern_grid = zeros(2, size(θ_grid)...)
    pattern::anten_pattern # for linear interpolation
    coeffi::ComplexF64  = 1
    #XXX maybe slow down the performance
    local_coord = Matrix(1.0I, 3, 3)
end

# # make it like a Vector{Float64}
# Base.getindex(x::anten_point, i::Int64) = x.p[i]
# Base.setindex!(x::anten_point, v::AbstractFloat, i) = x.p[i] = v
# Base.getindex(x::anten_point, I::Vararg{Int, N}) = Throw("cannot indexing over 1-d")
# Base.setindex!(x::anten_point, v::AbstractFloat, I::Vararg{Int, N}) = Throw("cannot indexing over 1-d")

function set_point_loc_coord!(point::anten_point, coord::Matrix{Float64})
    point.local_coord = coord
end

function set_point_pattern!(point::anten_point, pattern::anten_pattern)
    point.pattern = pattern
end

# generate points of N-element linear array
function point_linear(; N, dx, pattern)
    vec_2_struct_point = p -> anten_point(p=(x=p[1], y=p[2], z=p[3]), pattern=pattern)
    result = collect(zip(zeros(N), LinRange(0, (N - 1)dx, N), zeros(N)))
    vec_2_struct_point.(result)
end



function point_rectangle(; Nx, Ny, dx, dy, pattern)
    vec_2_struct_point = (p, coeffi) -> anten_point(p=p, pattern=pattern, coeffi=coeffi)
    point = [[dx * i, dy * j, 0] for i in 1:Nx for j in 1:Ny]

    # R0 = 10
    # n = 5
    # l = 3.5
    # I = Taylor_Chebyshev_I(n, R0)

    # coeffi = [I(l, x * 0.5) / I(l, 0) * I(l, y * 0.5) / I(l, 0)
    #           for x in (1:Nx) .- (Nx + 1) / 2
    #           for y in (1:Ny) .- (Ny + 1) / 2]
    # vec_2_struct_point.(point, coeffi)
    vec_2_struct_point.(point, 1)   
end


function create_array(; vec_points, pattern=pattern_identity)
    map(vec_points) do vec
        anten_point(p=vec, pattern=pattern)
    end
end




c = 299792458
# use Linrage instead of collect data can make interpolation.jl faster, collect it will make Interpolations.jl use non-uniform gird
θ_default::Vector{Float64}, ϕ_default::Vector{Float64} = (LinRange(0, pi, 180 +1), LinRange(-pi, pi, 360+1))
# const θ_default, ϕ_default = (LinRange(0, 180, 360 +1), [-180,  0, ]) .|> x -> deg2rad.(x)
θ_grid::Matrix{Float64}, ϕ_grid::Matrix{Float64} = ([θ for θ in θ_default, ϕ in ϕ_default],
    [ϕ for θ in θ_default, ϕ in ϕ_default])

# const cart_grid = [[sph2cart(θ,ϕ,1)...] for θ in θ_default, ϕ in ϕ_default]

# θ_default, ϕ_default = ([0,90], [0,0]) .|> x -> deg2rad.(x)


f::Float64 = 1
λ::Float64 = 1
k::Float64 = 1

function set_param(; freq)
    global f = freq
    global λ = c / f
    global k = 2pi / λ
    (f=f, λ=λ, k=k)
end

export
    anten_point,
    anten_pattern,
    θ_default, ϕ_default,
    θ_grid, ϕ_grid,
    # cart_grid,
    c,
    f,
    λ,
    k,
    set_param,
    set_point_loc_coord!,
    set_point_pattern!,
    create_array,
    pattern_identity,
    point_rectangle,
    pattern_dipole

