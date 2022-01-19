using CSV 
using antenna
using PlotlyJS
using Chain

# read file
raw_anten =filepath -> begin 
    CSV.File(filepath) |> DataFrame
end

A1 = raw_anten("test/A1.csv")
A2 = raw_anten("test/A2.csv")

phi = reshape(A1[:,1] .|> deg2rad, (181,91))
theta = reshape(A1[:,2].|> deg2rad, (181,91))

E_A1_phi = @chain A1[:,3] + A1[:,4]*1im reshape(_, (181,91))
E_A1_theta =@chain A1[:,5] + A1[:,6]*1im reshape(_, (181,91))
E_A2_phi = @chain A2[:,3] + A2[:,4]*1im reshape(_, (181,91))
E_A2_theta =@chain A2[:,5] + A2[:,6]*1im reshape(_, (181,91))

E_theta = E_A1_theta .+ E_A2_theta
E_phi = E_A1_phi .+ E_A2_phi
E = abs.(E_theta./1000).^2 .+ abs.(E_phi./1000).^2



tulple_matrix = sph2cart.(theta, phi, E)
x = [i.x for i = tulple_matrix]
y = [i.y for i = tulple_matrix]
z = [i.z for i = tulple_matrix]

@show max
plot(
    surface(x = x, y = y, z = z,),
    # Layout(
    #     scene_xaxis_range=[-max,max],
    #     scene_yaxis_range=[-max,max],
    #     scene_zaxis_range=[-max,max])
)



# _3Dvec_2_1Dvec = x -> LinRange(minimum(x), maximum(x), size(unique(x), 1))

# translate2Interpolation = (θ, ϕ, G) -> @chain begin
#     reshape(G, (size(ϕ, 1), size(θ, 1)))
#     permutedims((order[1], order[2]))
#     LinearInterpolation((θ, ϕ), _)
# end
# θ = _3Dvec_2_1Dvec(raw_anten.θ)
# ϕ = _3Dvec_2_1Dvec(raw_anten.ϕ)

# anten_pattern(
#     θ = translate2Interpolation(θ, ϕ, raw_anten.Gθ),
#     ϕ = translate2Interpolation(θ, ϕ, raw_anten.Gϕ)
# )

