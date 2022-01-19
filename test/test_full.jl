using BenchmarkTools
using Revise
using antenna
using LinearAlgebra

set_param(f = 2.6e9)

rotate45 = [1 0         0;
            0 1/sqrt(2) -1/sqrt(2);
            0 1/sqrt(2) 1/sqrt(2)]
            
rotate0 = Matrix(1.0I, (3, 3))

# test for wing
(vecx, vecy, vecz, point_wing) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")
vec_x_y_z_2_vec_coord_matrix = (vecx, vecy, vecz) -> [vecx vecy vecz]
coords = vec_x_y_z_2_vec_coord_matrix.(vecx, vecy, vecz)

print("\ec\n")
# read element file
elem_pattern = anten_read("test/pattern.csv")
# create point with pattern 
point = point_from_vec(vec_point = point_wing, pattern = elem_pattern)
# # set coord
set_point_loc_coord!.(point, coords)
# # calculate
θₜ, ϕₜ = (0, 90) .|> deg2rad
global_pattern =@time cal_pattern(point, θₜ, ϕₜ)
ploting_pattern(global_pattern,min=-20)
# ploting_point(point)



# # test for 8 linear
# point_wing = point_from_vec(vec_point = [[0,  0, 3e8/2.6e9/2*i] for i in 1:8], pattern=pattern_identity)
# ploting_point(point_wing)
# set_point_loc_coord!.(point_wing, [rotate45])
# θₜ, ϕₜ = (40, 90) .|> deg2rad
# global_pattern = @time cal_pattern(point_wing, θₜ, ϕₜ)
# # plot_E_vec_field(global_pattern)
# ploting_pattern(global_pattern,min=-10)
# # ploting_point(point_wing)



# # test for 8 x 8 planar array
# print("\ec\n")
# set_param(f=3e9)
# c = 299792458
# λ = c/3e9
# point_rec = point_rectangle(Nx=8, Ny=8,dx=λ/2, dy=λ/2, pattern = pattern_dipole)
# ploting_point(point_rec)
# θₜ, ϕₜ = (0, 80) .|> deg2rad
# global_pattern = cal_pattern(point_rec, θₜ, ϕₜ)
# ploting_pattern(global_pattern)