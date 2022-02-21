using BenchmarkTools
using Revise
using antenna
using LinearAlgebra
using PlotlyJS

set_param(f = 2.7e9)

print("\ec\n")
# read element file
elem_pattern = include("read_element_re.jl")
ploting_pattern(elem_pattern,min=-35)
radiated_power(elem_pattern)
directivity(elem_pattern) |> maximum |> x-> 10log10(x)

radiation_intensity(elem_pattern) |> maximum
gain(elem_pattern, 0.948039) |> maximum |> x-> 10log10(x)

