using Antenna
using Meshes
using Rotations
using GLMakie
using LinearAlgebra

function facet_filter!(vec_facet)
    filter!(vec_facet) do x
        begin 
            y_lim = -3.7
            z_lim = 1.385
            # cut y
            x.vertex[1][2] < y_lim && x.vertex[2][2] < y_lim &&x.vertex[3][2] < y_lim &&
            # cut z
            x.vertex[1][3] < z_lim && x.vertex[1][3] < z_lim && x.vertex[1][3] < z_lim
        end 
    end
end
function create_plane_mesh(vec_facet)
    rotate_facet(vec_facet)
    facet_filter!(vec_facet)
    plane_mesh = vecfacet2ply(vec_facet)
end



set_param(freq=2.7e9)
## ---------test
patch = read_hfss_pattern("test/elem/rE.csv")
vec_facet = create_facets_from_file()
plane_mesh = create_plane_mesh(vec_facet)

## sample
plane_sample = sample(plane_mesh, MinDistanceSampling(λ*0.7*1.01, ρ=0.95)) |>collect 
array = create_array_from_meshes_sample(plane_sample, plane_mesh, pattern=patch)
scatter(array)
plot_pattern(res)
fig = Figure(resolution=(1920,1080))
display(fig)
ax, sc = scatter(fig[1,1],array, color = abs.(D))
res,D = cal_pattern(array, 0,0, spin=true, optimal_directivity=:all);
ax, sc_pattern = pattern3d(fig[1,2], res)
slider1 = Slider(fig[2,1], range=range(0,-pi,100),startvalue=-pi/2)
# angle_scan = mod_angle_rad.(pi/2, range(-pi,0,100)) |>stack

capture_pattern = nothing
lift(slider1.value) do i
    θ,ϕ = mod_angle_rad.(pi/2, i) 
    res,D = cal_pattern(array, θ,ϕ, spin=true, optimal_directivity=:all);
    # println(directivity_beta(res, θ,ϕ) |> x->10log10(x)) 
    global capture_pattern = res
    println(directivity_beta(res, θ,ϕ) )
    sc.color[] = abs.(D)
    sc_pattern[1][] = res
    sleep(0.02)
    @show θ,ϕ
end
res
res = directivity_beta(res, θ_grid, ϕ_grid)
pattern3d(res)
##
# ##
# fig = Figure(viewmode=:fit, aspect=(1,1,1))
# display(fig)
# ax = Axis3(fig[1,1], aspect=(1,1,1))
# slider1 = Slider(fig[2,1], range=range(0.5,1,100),startvalue=0.6)
# res = sample(plane_mesh, MinDistanceSampling(λ*0.5*1.01, ρ=0.6)) |>collect 
# t = Observable(res)
# sc = scatter!(ax, t) 

# lift(slider1.value) do i
#     @show i
#     @time res = sample(plane_mesh, MinDistanceSampling(λ*i*1.00, ρ=0.9)) |>collect 
#     t[] = res
#     t[] = deleteat!(t[], length(t.val))
#     t[] = push!(t[], [res[1].coords...])   
#     sleep(0)
# end
