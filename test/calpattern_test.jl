using Revise
using Antenna
using Plots
using BenchmarkTools
using Tullio

dipole = read_hfss_pattern()
set_param(freq=2.7e9)
point = point_rectangle(Nx = 10,Ny=10, dx=λ/2, dy=λ/2,pattern=dipole)
@time res_pattern = cal_pattern(point, 0.,0., spin=false);
plot(res_pattern)



