using DataFrames
using CSV
using Chain
using Antenna
using Interpolations

print("\ec")
extract = (filename) -> @chain CSV.File(filename)  DataFrame _[181:end, 2:end] Matrix circshift(_, (0, -180))
extract_theta_phi = (filename) -> @chain CSV.File(filename)  DataFrame _[181:end, 2:end] Matrix circshift(_, (0, -180))
# theta = @chain CSV.File("test/elem/mag_theta.csv")  DataFrame _[181:end,1] deg2rad.(_)
theta = collect(0.:180.) .|> deg2rad
phi = collect(-180.:180.) .|> deg2rad 

mag_theta = extract("test/elem/mag_theta.csv")
mag_phi = extract("test/elem/mag_phi.csv")
ang_theta = extract("test/elem/ang_theta.csv")
ang_phi = extract("test/elem/ang_phi.csv") 

rE_theta = [mag*exp(-1im*deg2rad(ang)) for (mag,ang) in zip(mag_theta, ang_theta)]
rE_phi = [mag*exp(-1im*deg2rad(ang)) for (mag,ang) in zip(mag_phi, ang_phi)]

anten_pattern(
    θ=LinearInterpolation((theta, phi), rE_theta),
    ϕ=LinearInterpolation((theta, phi), rE_phi)
)