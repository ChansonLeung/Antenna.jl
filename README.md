# Antenna.jl

# Introduction


Antenna.jl is a julia package that helps you design the antenna array. It's a personal project that aims to facilitate prototyping the custom type of antenna array like conformal array and better studying the properties and optimization about antenna array. It has the below features

1. Synthesis. You can calculate the performance of the antenna array with imported element antenna data. You can either use the data from simulation software or test instrument.  When dealing with the conformal array, the polarization, element rotation is in consideration.
2. Plotting. It is based on the PlotlyJS. you can plot the directivity pattern, gain pattern, in both 3D and 2D modes, and calculate some parameters like radiated power, radiated efficiency, etc.
3. Light weight short code.

This project is still in the initial stage. Now the main goals are

1. Optimization, which can be able to create your optimization model and bind it to the designed antenna array.
2. Acceptable performance, which makes use of the parallel properties of the julia language to keep all the calculation fast enough without involving the vectorization code like Matlab or Numpy.
3. Friendly usage with good documentation.

# Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ChansonLeung/Antenna.jl")
```

# Example
