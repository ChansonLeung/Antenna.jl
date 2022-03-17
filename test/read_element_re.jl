using DataFrames
using CSV
using Chain
using Interpolations
using Antenna

extract_theta_phi = (filename) -> @chain begin
    CSV.File(filename)  
    DataFrame 
    Matrix 
end
data_csv = extract_theta_phi("test/elem/rE.csv")
theta = reshape(data_csv[:,1], (361,181))
phi = reshape(data_csv[:,2], (361,181))

rE_phi = data_csv[:, 3] + data_csv[:, 4]*1im |>
    (x -> reshape(x, (361,181))) .|>
    (x -> x/1000)  |>
    x -> permutedims(x, [2,1])

rE_theta = data_csv[:, 5] + data_csv[:, 6]*1im |> 
    x -> reshape(x, (361,181)) .|>
    (x -> x/1000) |>
    x -> permutedims(x, [2,1])




theta = collect(0.:180.) .|> deg2rad
phi = collect(-180.:180.) .|> deg2rad 


return anten_pattern(
    θ=LinearInterpolation((theta, phi), rE_theta),
    ϕ=LinearInterpolation((theta, phi), rE_phi)
)