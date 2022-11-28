using Revise
using Antenna

elem_pattern = include("read_element_re.jl")

plot_pattern(elem_pattern,min=-35)


set_param(f=2.7e9)

point = point_rectangle(Nx = 8,Ny=8, dx=λ/2, dy=λ/2,pattern=pattern_identity)
plot_point(point)

res_pattern = cal_pattern(point, 0,0)

# @show radiated_power(res_pattern)

plot_pattern(res_pattern, min=-30)

plot_pattern_2D(res_pattern) 