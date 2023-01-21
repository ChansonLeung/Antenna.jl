using CUDA
using LinearAlgebra
using BenchmarkTools
using Tullio
using StaticArrays
using FLoops

function sph2cart(θ, ϕ, r)
    @SVector[r * sin(θ)cos(ϕ),
        r * sin(θ)sin(ϕ),
        r * cos(θ)]
end

function cart2sph(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    @SVector [acos(z / r),
        atan(y, x),
        r]
end

function rotation_kernel(grid_in, grid_out, θ_grid, ϕ_grid, mat_rot)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for i in idx:stride:length(length(θ_grid))
        # rotate grid
        θ, ϕ = θ_grid[i], ϕ_grid[i]
        sinθ, sinϕ, cosθ, cosϕ = sin(θ), sin(ϕ), cos(θ), cos(ϕ)
        # FIXME need matmul
        x, y, z = 1 * @SVector[sinθ * cosϕ, sinθ * sinϕ, cosθ]
        θ_rot, ϕ_rot, r_rot = cart2sph(x, y, z)

        # get pattern vector value in rotate grid
        sinθ_rot, sinϕ_rot, cosθ_rot, cosϕ_rot = sin(θ_rot), sin(ϕ_rot), cos(θ_rot), cos(ϕ_rot)
        # FIXMEE need interpolation and matmul
        val_vec_θ_rot = 1.0 * mat_rot' * @SVector[cosθ_rot * cosϕ_rot, cosθ_rot * sinϕ_rot, -sinθ_rot]
        val_vec_ϕ_rot = 1.0 * mat_rot' * @SVector[-sin(ϕ_rot), cos(ϕ_rot), 0.0]

        # project vector to global grid
        vec_θ_global = @SVector[cosθ * cosϕ, cosθ * sinϕ, -sinθ]
        vec_ϕ_global = @SVector[-sinϕ, cosϕ, 0.0]

        # FIXME Why out of bounds?
        grid_out[1, (i-1)%361+1, (i-1)÷361+1, blockIdx().y] = dot(val_vec_θ_rot, vec_θ_global) + dot(val_vec_ϕ_rot, vec_θ_global)
        grid_out[2, (i-1)%361+1, (i-1)÷361+1, blockIdx().y] = dot(val_vec_ϕ_rot, vec_ϕ_global) + dot(val_vec_θ_rot, vec_ϕ_global)
    end
    return nothing
end

const θ = LinRange(0, pi, 361)
const ϕ = LinRange(-pi, pi, 181)
const θ_grid = [θ for θ in θ, ϕ in ϕ]
const ϕ_grid = [ϕ for θ in θ, ϕ in ϕ]
rot_size = (361, 181)
points_num = 100
grid_in = rand(2, rot_size..., points_num)
grid_out = rand(2, rot_size..., points_num)
mat_rot = rand(3, 3)

@time A =  rand(2,361,181,100);
@time Array(cu(A));

f(grid_in, grid_out, θ_grid, ϕ_grid, mat_rot) = begin
    # allocate GPU array
    begin
        grid_in = cu(grid_in)
        grid_out = cu(grid_out)
        θ_grid = cu(θ_grid)
        ϕ_grid = cu(ϕ_grid)
        mat_rot = cu(SMatrix{3,3}(mat_rot))
    end
    # launch kernel
    threads_num = 128
    blocks_num = (ceil(Int, length(grid_in) / threads_num), points_num)
    CUDA.@sync @cuda threads = threads_num blocks = blocks_num rotation_kernel(grid_in, grid_out, θ_grid, ϕ_grid, mat_rot)
    Array(grid_out)
end

@time f(grid_in, grid_out, θ_grid, ϕ_grid, mat_rot);

