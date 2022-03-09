using BenchmarkTools
using Revise
using antenna

set_param(f = 2.7e9)

# test for wing
# (vecx, vecy, vecz, point_wing) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")
# vec_x_y_z_2_vec_coord_matrix = (vecx, vecy, vecz) -> [vecx vecy vecz]
# coords = vec_x_y_z_2_vec_coord_matrix.(vecx, vecy, vecz)

# (_,_,_, point_wing_0) = include("D:/ChansonDocument/project/2112_7_机翼工况计算/yang/read_data.jl")
# point_0 = point_from_vec(vec_point = point_wing_0, pattern = elem_pattern) 


print("\ec\n")
# read element file
elem_pattern = include("read_element_re.jl")
ploting_pattern(elem_pattern,min=-35)

# create point with pattern 
# point = point_from_vec(vec_point = point_wing, pattern = elem_pattern)
# ploting_point(point)
# # set coord
set_point_loc_coord!.(point, coords)
# # calculate
# θₜ, ϕₜ = (0, 90) .|> deg2rad
# global_pattern =@time cal_pattern(point,point_0, θₜ, ϕₜ)

# fig = ploting_pattern(global_pattern)
# ploting_pattern(elem_pattern)
# file = open("0.html", "w")
# PlotlyJS.savefig(file, fig,format = "html")
# close(file)
# fig = ploting_point(point)

