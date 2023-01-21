using Revise
using Antenna
using Plots

@time dipole = read_hfss_pattern();
set_param(freq=2.7e9)
@time point = point_rectangle(Nx = 10,Ny=10, dx=λ/2, dy=λ/2,pattern=dipole);
res_pattern = cal_pattern(point, 0.,0.);
plot(res_pattern)
