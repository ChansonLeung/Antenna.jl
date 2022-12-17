using Antenna
using DataFrames
using CSV
using Interpolations

export 
    read_hfss_pattern

function read_hfss_pattern(name = "test/dataset/rE_theta_phi_dipole.csv", mode=:polar)

    data_csv = CSV.read(name, DataFrame)

    if mode == :polar
        # 1     2      3            4           5           6
        # phi theta phi_re[mv] phi_im[mv] theta_re[mv] theta_im[mv]

        (phi_grid,
        theta_grid,
        rE_phi_re,
        rE_phi_im,
        rE_theta_re,
        rE_theta_im) = getindex.([data_csv], !,["Phi[deg]","Theta[deg]","re(rEPhi)[mV]","im(rEPhi)[mV]","re(rETheta)[mV]","im(rETheta)[mV]"]) .|> 
                    x->Vector(Float64.(x)) |>
                    x->reshape(x,361,181)


        rE_phi = rE_phi_re + rE_phi_im * 1im |>
                x -> permutedims(x, [2, 1])/1000

        rE_theta = rE_theta_re + rE_theta_im * 1im |> 
                x -> permutedims(x, [2, 1])/1000


        theta = collect(0.0:180.0) .|> deg2rad 
        phi = collect(-180.0:180.0) .|> deg2rad 



        return anten_pattern(
            θ=linear_interpolation((theta, phi), rE_theta, extrapolation_bc=Flat()),
            ϕ=LinearInterpolation((theta, phi), rE_phi, extrapolation_bc=Flat())
        )
    else mode == :xyz
        #   1           2               3               4               5               6               7               8
        # Phi[deg], Theta[deg],     re(rEX)[mV],    im(rEX)[mV],    re(rEY)[mV],    im(rEY)[mV],    re(rEZ)[mV],    im(rEZ)[mV]
        data_csv = extract_theta_phi(name) |>
                x -> reshape(x, 361, 181, :)
        (
            (phi_grid,
            theta_grid,
            rE_ex_re,
            rE_ex_im,
            rE_ey_re,
            rE_ey_im,
            rE_ez_re,
            rE_ez_im) = collect(eachslice(data_csv, dims=3))
        )

        rE_ex = rE_ex_re + rE_ex_im * 1im |>
                x -> permutedims(x, [2, 1])/1000

        rE_ey = rE_ey_re + rE_ey_im * 1im |>
                x -> permutedims(x, [2, 1])/1000
        rE_ez = rE_ez_re + rE_ez_im * 1im |>
                x -> permutedims(x, [2, 1])/1000


    end
end

