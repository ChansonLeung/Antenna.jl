# FIXME
module type
using LinearAlgebra
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
pattern_identity = anten_pattern(θ = (θ, ϕ) -> 1, ϕ = (θ, ϕ) -> 0)
# antenna theory 4-84
pattern_dipole = anten_pattern(
    θ = (θ, ϕ) -> 1im*η*exp(-1im*k)/(2pi) * (cos(pi/2*cos(θ))/(sin(θ)+1e-6)),
    ϕ = (θ, ϕ) -> 0
)


Base.@kwdef mutable struct anten_point
    #point
    p::NamedTuple{(:x, :y, :z),Tuple{Float64,Float64,Float64}}
    pattern::anten_pattern
    #XXX maybe slow down the performance
    local_coord = Matrix(1.0I, 3, 3)
end

function set_point_loc_coord!(point::anten_point, coord::Matrix{Float64})
    point.local_coord = coord
end
function set_point_pattern!(point::anten_point, pattern::anten_pattern)
    point.pattern = pattern
end

# generate points of N-element linear array
function point_linear(; N, dx, pattern)
    vec_2_struct_point = p -> anten_point(p = (x = p[1], y = p[2], z = p[3]), pattern = pattern)
    result = collect(zip(zeros(N), LinRange(0, (N - 1)dx, N), zeros(N)))
    vec_2_struct_point.(result)
end

function point_rectangle(;Nx,Ny, dx,dy,pattern)
    vec_2_struct_point = p -> anten_point(p = (x = p[1], y = p[2], z = p[3]), pattern = pattern)
    point = [[dx*i, dy*j, 0] for i in 1:Nx for j in 1:Ny]
    vec_2_struct_point.(point)   
end

function point_from_vec(; vec_point, pattern)
    vec_2_struct_point = (vec) -> anten_point(p = (x = vec[1], y = vec[2], z = vec[3]), pattern = pattern)
    vec_2_struct_point.(vec_point)
end




c = 299792458
const θ_default, ϕ_default = (LinRange(0,180,361), LinRange(-180,180,361)) .|> x -> deg2rad.(x)
const θ_grid, ϕ_grid = ([θ for θ in θ_default, ϕ in ϕ_default],
                 [ϕ for θ in θ_default, ϕ in ϕ_default])

# θ_default, ϕ_default = ([0,90], [0,0]) .|> x -> deg2rad.(x)


f = 1
λ = 1
k = 1

function set_param(; f)
    λ = c / f
    global k = 2pi / λ
    (f = f, λ = λ, k = k)
end
function set_grid(θ,ϕ)
    global θ_grid = θ
    global ϕ_grid = ϕ
end
export
    anten_point,
    anten_pattern,
    θ_default, ϕ_default,
    θ_grid, ϕ_grid,
    c,
    f,
    λ,
    k,
    set_param,
    set_grid,
    set_point_loc_coord!,
    set_point_pattern!,
    point_from_vec,
    pattern_identity,
    point_rectangle,
    pattern_dipole

end
