U₀ = 1/4pi
η = 120pi
E = (θ,ϕ) -> 1im*η*exp(-1im*k)/(2pi) * (cos(pi/2*cos(θ))/(sin(θ)+1e-6))
W = (θ,ϕ) -> 1/2η * abs(E(θ,ϕ))^2
θ = LinRange(0,pi, 100)

[W(θ,ϕ) for θ = θ, ϕ = 0]