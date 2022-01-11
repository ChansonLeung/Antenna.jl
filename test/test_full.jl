using BenchmarkTools
using Revise
using antenna
using LinearAlgebra

(vecx, vecy, vecz, point_wing) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")
set_param(f = 2.6e9)

# print("\ec")
print("\ec")
elem_pattern = anten_read("test/pattern.csv")
point = point_from_vec(vec_point = point_wing, pattern = elem_pattern)
vec_x_y_z_2_vec_coord_matrix = (vecx, vecy, vecz) -> [vecx vecy vecz]

coords = vec_x_y_z_2_vec_coord_matrix.(vecx, vecy, vecz)
rotate45 = [1 0         0;
            0 1/sqrt(2) -1/sqrt(2);
            0 1/sqrt(2) 1/sqrt(2)]
            
rotate0 = Matrix(1.0I, (3, 3))

set_point_loc_coord!.(point, coords)

θₜ, ϕₜ = (90, 90) .|> deg2rad

# global_pattern = @time cal_pattern(point, θₜ, ϕₜ)
# ploting_point(point)

point_wing = point_from_vec(vec_point = [[0,  0,3e8/2.6e9/1.1*i] for i in 1:8], pattern=pattern_identity)
set_point_loc_coord!.(point_wing, [rotate0])
ploting_point(point_wing)
global_pattern = @time cal_pattern(point_wing, θₜ, ϕₜ)
# global_pattern = rotate_pattern(pattern_identity, rotate45)

# plot_E_vec_field(global_pattern)
ploting_pattern(global_pattern)


# plot(
#     cone(
#         x = x,
#         y = y,
#         z = z,
#         u = u,
#         v = v,
#         w = w,
#         sizemode = "scaled",
#         sizeref = 30
#     )
# )
# A = 2:15
# B = 3:10
# [A for 45 in B for A in A ] ==vec([A for A in A , B in B])
