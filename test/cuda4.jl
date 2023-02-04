## 
using CUDA
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using Antenna
using Plots
using Tullio

## 
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
AF(point, θ, ϕ, k) = exp(1im * k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k)^-1

function cal_pattern_cuda_kernel(
    positions, 
    coeffi, 
    mat_rot, 
    pattern_in_re, 
    pattern_in_im, 
    pattern_out, 
    k,
    θₜ=0, 
    ϕₜ=0, 
    θ_default=cu(LinRange(0, pi, 181)), 
    ϕ_default=cu(LinRange(-pi, pi, 361)),
)
    idx_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride_i = gridDim().x * blockDim().x

    for i in idx_i:stride_i:length(θ_default)*length(ϕ_default)
        for j in 1:length(coeffi)
            mat_rot_i = SMatrix{3,3}(mat_rot[j])
            idx_θ = (i - 1) % length(θ_default) + 1
            idx_ϕ = (i - 1) ÷ length(θ_default) + 1

            # rotate grid
            θ, ϕ = θ_default[idx_θ], ϕ_default[idx_ϕ]
            sinθ, sinϕ, cosθ, cosϕ = sin(θ), sin(ϕ), cos(θ), cos(ϕ)
            x, y, z = mat_rot_i' * @SVector[sinθ * cosϕ, sinθ * sinϕ, cosθ]
            θ_rot, ϕ_rot, r_rot = cart2sph(x, y, z)

            # get pattern vector value in rotate grid
            sinθ_rot, sinϕ_rot, cosθ_rot, cosϕ_rot = sin(θ_rot), sin(ϕ_rot), cos(θ_rot), cos(ϕ_rot)
                # FIXME need interpolation
            tex_idx_θ = θ_rot * (length(θ_default)-1)/pi+1
            tex_idx_ϕ = (ϕ_rot +pi) * (length(ϕ_default)-1)/2pi+1
            Gθ_rot_local = pattern_in_re[1,tex_idx_θ, tex_idx_ϕ] + 1im*pattern_in_im[1,tex_idx_θ, tex_idx_ϕ]
            Gϕ_rot_local = pattern_in_re[2,tex_idx_θ, tex_idx_ϕ] + 1im*pattern_in_im[2,tex_idx_θ, tex_idx_ϕ]

            val_vec_θ_rot =  Gθ_rot_local* mat_rot_i * @SVector[cosθ_rot * cosϕ_rot, cosθ_rot * sinϕ_rot, -sinθ_rot]
            val_vec_ϕ_rot =  Gϕ_rot_local* mat_rot_i * @SVector[-sinϕ_rot, cosϕ_rot, 0.0]

            # project vector to global grid
            vec_θ_global = @SVector[cosθ * cosϕ, cosθ * sinϕ, -sinθ]
            vec_ϕ_global = @SVector[-sinϕ, cosϕ, 0.0]

            G_θ = dot(val_vec_θ_rot, vec_θ_global) + dot(val_vec_ϕ_rot, vec_θ_global)
            G_ϕ = dot(val_vec_ϕ_rot, vec_ϕ_global) + dot(val_vec_θ_rot, vec_ϕ_global)

            sinθₜ, sinϕₜ, cosθₜ, cosϕₜ = sin(θₜ), sin(ϕₜ), cos(θₜ), cos(ϕₜ)
            p = @SVector[positions[1, j], positions[2, j], positions[3, j]]
            # # ---
            # tmp = coeffi[j] * Iₛ(p, θₜ, ϕₜ, k) * AF(p, θ, ϕ, k)
            tmp = coeffi[j] * exp(1.0im * k * (p[1] * sinθ*cosϕ + p[2] * sinθ*sinϕ + p[3] * cosθ))
            # exp(1.0im * k * (p[1] * sinθₜ*cosϕₜ + p[2] * sinθₜ*sinϕₜ + p[3] * cosθₜ))

            pattern_out[1, idx_θ, idx_ϕ] += tmp * G_θ
            pattern_out[2, idx_θ, idx_ϕ] += tmp * G_ϕ
        end
    end
    return nothing
end

##
set_param(freq=2.7e9)
pattern = read_hfss_pattern()
pattern_size = (181, 361)
points_num = 100


using Rotations
mat_rot = fill(SMatrix{3,3}((1.0f0*I(3)*RotX(pi/4)*RotZ(pi/4) )), points_num)
positions = begin
    [@SVector[λ/2*i, λ/2*j, 0] for i in 1:10 for j in 1:10]
end |> x->reinterpret(reshape, Float64, x)
coeffi = ones(ComplexF64,points_num).*Iₛ.(eachcol(positions), 0, 0, k)
pattern_out_cu = CUDA.zeros(ComplexF32, 2,pattern_size...)

pattern_in = [
    pattern.θ.(θ_grid,ϕ_grid) |> x->reshape(x, 1,size(x)...);
    pattern.ϕ.(θ_grid,ϕ_grid) |> x->reshape(x, 1,size(x)...);
]
heatmap(pattern_in[1,:,:].|>abs)

pattern_in_re = CuTextureArray(cu(real.(pattern_in))) |> x -> CuTexture(x; interpolation=CUDA.LinearInterpolation())
pattern_in_im = CuTextureArray(cu(imag.(pattern_in))) |> x -> CuTexture(x; interpolation=CUDA.LinearInterpolation())

f() = @cuda threads = 128 blocks=ceil(Int, 361*181/128) cal_pattern_cuda_kernel(
    cu(positions),
    cu(coeffi),
    cu(mat_rot),
    pattern_in_re,
    pattern_in_im,
    pattern_out_cu,
    Antenna.k
);
f()
pattern_out = Array(pattern_out_cu)

res = abs.(pattern_out)
res = @. res[1,:,:]^2 + res[2,:,:]^2
heatmap(res)

## calculate conformal, 随便设的随机点
antenna_array = begin
    local p = Vector{anten_point}()
    for (p_i, r_i) in zip(eachslice(positions,dims=2), mat_rot)
        array_element = anten_point(p=p_i, pattern=pattern, local_coord=r_i);
        push!(p, array_element)
    end
    p
end
result_pattern = cal_pattern(antenna_array, 0.0, 0.0, spin=true);
plot(result_pattern)

begin
    Gθ =result_pattern.θ(θ_default,ϕ_default) .|>abs
    Gϕ =result_pattern.ϕ(θ_default,ϕ_default) .|>abs
    @. Gθ^2 + Gϕ^2
# end .- res
end ≈ res
# end |>heatmap

