using Antenna
using Plots
using Rotations
using FLoops
using BenchmarkTools

dipole_hfss = (
    read_hfss_pattern("test/dataset_accurate/rE_theta_phi_dipole.csv"),
    read_hfss_pattern("test/dataset/rE_theta_phi_dipole.csv")
)
dipole_hfss_rot = (
    read_hfss_pattern("test/dataset_accurate/rE_theta_phi_dipole_rotate.csv"),
    read_hfss_pattern("test/dataset/rE_theta_phi_dipole_rotate.csv")
)
dipole_rot = (
    rotate_pattern_tullion(zeros(ComplexF64, 2, 181, 361), dipole_hfss[1], RotZ(pi / 4) * RotY(pi / 4) |> Matrix),
    rotate_pattern_tullion(zeros(ComplexF64, 2, 181, 361), dipole_hfss[2], RotZ(pi / 4) * RotY(pi / 4) |> Matrix)
);

function pattern_error_abs(pattern1, pattern2, field)
    G_θ1 = getfield(pattern1, field).(θ_grid, ϕ_grid)
    G_θ2 = getfield(pattern2, field).(θ_grid, ϕ_grid)

    (@. (abs(G_θ1) - abs(G_θ2)),
    @. (abs(G_θ1) - abs(G_θ2)) / abs(G_θ2))
end

(error_abs, error_rel) = pattern_error_abs(dipole_hfss_rot[2], dipole_rot[2], :θ)
(error_abs_accr, error_rel_accr) = pattern_error_abs(dipole_hfss_rot[1], dipole_rot[1], :θ)

plot(
    histogram(vec(error_abs), title="幅度绝对误差(HFSS 2°采样)", color=:red),
    histogram(vec(error_rel), title="幅度相对误差(HFSS 2°采样)", xlims=(-0.1, 0.1), color=:red),
    histogram(vec(error_abs_accr), title="幅度绝对误差(HFSS 1°采样)"),
    histogram(vec(error_rel_accr), title="幅度相对误差(HFSS 1°采样)", xlims=(-0.1, 0.1)),
    legend=false,
    size=(800, 600)
)
plot(
    plot(dipole_hfss_rot[1], title="HFSS建模旋转"),
    plot(dipole_rot[1], title="自定义旋转"),
    plot(dipole_hfss[1], title="HFSS原始"),
    size=(800, 600),
    layout=(2, 2)
)
