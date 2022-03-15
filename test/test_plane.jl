using Revise

include("../src/antenna.jl")
using .antenna


elem_pattern = include("read_element_re.jl")
ploting_pattern(elem_pattern,min=-35)

set_param(f=2.7e9)
point = point_rectangle(Nx = 8,Ny=8, dx=0.111/2, dy=0.111/2,pattern=elem_pattern)
res_pattern = cal_pattern(point, point, 0,0)
ploting_pattern(res_pattern, min=-20)
# ploting_point(point)