AF(point, θ, ϕ, k) = exp(1im * k * (point[1] * sin(θ)cos(ϕ) + point[2] * sin(θ)sin(ϕ) + point[3] * cos(θ)))
Iₛ(point, θₜ, ϕₜ, k) = AF(point, θₜ, ϕₜ, k) * exp(-1)
phase_move2coplane = (point::Vector{Float64},f::Float64, θₜ::Float64, ϕₜ::Float64) ->begin
    c = 299792458
    k = 2pi*f/c
    AF(point, θₜ, ϕₜ, k)*Iₛ(point, θₜ, ϕₜ, k) |> angle
end

# I_compensate 计算阵列天线中一个点弯曲前后需要的补偿相位
# point_bent: 机翼弯曲后的天线坐标点位置 in {x,y,z} m
# I_flat:     机翼未形变时的馈电相位 in Rad
# f:          工作频率 in Hz
# θₜ, ϕₜ:       天线的扫描角度 in Rad
# return -> 每个单元需要补偿的相位量 in Rad
function I_compensate(point_bent::Vector{Float64} , I_flat::Float64, f::Float64, θₜ::Float64, ϕₜ::Float64)
    I_bent = phase_move2coplane(point_bent, f, θₜ, ϕₜ)
    I_compensate = I_bent.- I_flat
end

# 对所有点进行补偿
function I_compensate_all(point_bent::Vector{Vector{Float64}}, 
                            point_phase_origin::Vector{Float64}, 
                            I_flat::Vector{Float64},
                            f::Float64, 
                            θₜ::Float64, 
                            ϕₜ::Float64)
    I_compensate_all =  (I_compensate.(point_bent, I_flat, f, θₜ, ϕₜ) .-
                        phase_move2coplane(point_phase_origin,f,θₜ,ϕₜ)) .|> mod2pi
end


point_bend, point_phase_origin, I_flat = include("read_data.jl")
f=2.4
θ = 0.
ϕ = 0.

using BenchmarkTools
@btime I_compensate_all(point_bend, point_phase_origin, I_flat, f,θ,ϕ)



