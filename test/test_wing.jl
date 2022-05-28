using Revise
using Antenna
using BenchmarkTools
using FileIO

target = (0.0, 0.0)

# read wing shape data
file_wing = load("test/local_coordinate_and_point.jld2")
points, local_coordinate = file_wing["-400"]
points_0, local_coordinate_0 = file_wing["0"]

plot_point(points)

# read antenna data
elem_pattern = include("read_element_re.jl")
# set frequency
set_param(f=2.7e9)

# set point and coordinate
array = create_array(vec_points=points_0, pattern=elem_pattern)
# calculate result
res_pattern = cal_pattern(array, target...)

# plot result
SLL , max_gain = SLL_maxgain(res_pattern,ϕ=90) .|> x->10log10(x)
println("副瓣电平：$(SLL)dB, 最大方向性系数: $(max_gain)dB ")
plot_pattern_all(res_pattern,θ=-90:0.1:90, ϕ = 90, point = points_0)
