using Revise
using Antenna
using BenchmarkTools
using Plots
target = (0.0,0.0)
set_param(f=24e9)
c = 299792458
f = 24e9
λ = c / f

array = point_rectangle(Nx=10, Ny=10,dx=λ/2, dy=λ/2, pattern = pattern_identity)
res =@time cal_pattern(array, target...)
@time SLL_maxgain(res)
# plot_pattern(res)

