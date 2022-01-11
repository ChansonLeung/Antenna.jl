using Plotly
using BenchmarkTools
using LinearAlgebra
using Revise
using antenna

(vecx,vecy,vecz,point_wing) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")

print("\ec")

set_param(f=2.6e9)

elem_pattern = anten_read("test/pattern.csv")

# point = point_linear(N = 8, dx = 25e-3, pattern = elem_pattern)

point = point_from_vec(vec_point = point_wing, pattern=elem_pattern)
# point = point_from_vec(vec_point = point_wing, pattern=pattern_identity)

vec_x_y_z_2_vec_coord_matrix = (vecx,vecy,vecz) -> [vecx vecy vecz]
coords = vec_x_y_z_2_vec_coord_matrix.(vecx,vecy,vecz)
# set_point_loc_coord!.(point, [Matrix(1.0I, (3, 3))])
set_point_loc_coord!.(point, coords)

θₜ, ϕₜ = (0, 90) .|> deg2rad

# global_pattern = cal_pattern(point, θₜ, ϕₜ)
global_pattern =@time cal_pattern(point, θₜ, ϕₜ) # pattern

# heatmap([sqrt(global_pattern.θ(θ,ϕ)^2 + global_pattern.ϕ(θ,ϕ)^2) for θ=antenna.θ_default, ϕ=antenna.ϕ_default])




# maybe the problem of component combining