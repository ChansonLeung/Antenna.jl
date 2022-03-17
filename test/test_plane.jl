using Revise
using Antenna

elem_pattern = include("read_element_re.jl")
ploting_pattern(elem_pattern,min=-35)
set_param(f=2.7e9)
point = point_rectangle(Nx = 8,Ny=8, dx=0.111/2, dy=0.111/2,pattern=elem_pattern)
# point = point_rectangle(Nx = 4,Ny=4, dx=0.111/2, dy=0.111/2,pattern=elem_pattern)
res_pattern = cal_pattern(point, point, 0,0)
@show radiated_power(res_pattern)
ploting_pattern(res_pattern, min=-30)
Antenna.ploting_pattern_2D(res_pattern)