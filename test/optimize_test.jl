using Antenna
using Makie
using GLMakie
using Observables


patch = read_hfss_pattern("test/elem/rE.csv")
set_param(freq=2.7e9)


function mk_expand_points(point)
    x,y,z = [getindex.(point, i) for i in 1:3]
end

anten_get_points(array) = begin
    local x = getfield.(array, :p) .|> x->x[1]
    local y = getfield.(array, :p) .|> x->x[2]
    local z = getfield.(array, :p) .|> x->x[3]

    local coords = getfield.(array, :local_coord)
    local u = map(x->x[1,1], coords)
    local v = map(x->x[2,1], coords)
    local w = map(x->x[3,1], coords)
    [x,y,z, u, v, w]

end 
## ------------------calculate-----------------------------
 res_pattern, D = @time cal_pattern(array_sphere, deg2rad(0.), 0.0, spin=true, optimal_directivity=:all);
#  cal_pattern(array_sphere, deg2rad(0.), 0.0, spin=true, optimal_directivity=:all);
plot_pattern(res_pattern)



## --------------------plot------------------------------

fig = Figure()
ax = Axis3(fig[1,1])

slider_θ = Slider(fig[2,1], range=-90:90, startvalue=0)
slider_ϕ = Slider(fig[3,1], range=-90:90, startvalue=0)
begin
    local x,y,z =  mk_expand_points(p_sphere)
    # surface(x,y,z, color=ones(length(x)))
    surface!(x,y,z, color=fill(:gray, size(z)...))
end

color_points =  Observable(abs.(D))
begin 
    local x,y,z,o = anten_get_points(array_sphere)
    scatter!(x,y,z, color=color_points)
end

lift(slider_θ.value, slider_ϕ.value) do θ,ϕ
    θ,ϕ = mod_angle_rad(deg2rad(θ), deg2rad(ϕ))
    @time res_pattern, D = cal_pattern(array_sphere, θ,ϕ, spin=true, optimal_directivity=:all);
    color_points[] = abs.(D).^2
    sleep(0)
end
fig

## ---------------plot pattern in real time-------------------------
fig = Figure(resolution=(2000,1000))
ax_points = Axis3(fig[1,2], viewmode=:fit, aspect=:data)
ax = Axis3(fig[1,1], viewmode=:fit, aspect=:data)

slider_θ = Slider(fig[2,:], range=-180:180, startvalue=0)
slider_ϕ = Slider(fig[3,:], range=-180:180, startvalue=0)
x,y,z,r = [Observable(zeros(size(θ_grid))) for _ in 1:4]
D = Observable(zeros(length(array_sphere)))
surface!(ax,x,y,z,color=r, colormap="lightrainbow")
surface!(ax,x,y,z,color=r)
begin
    local x,y,z = [getindex.(p, i) for i in 1:3]
    surface!(ax_points, x,y,z, color=fill(:gray, size(x)...))
    local x,y,z =  getfield.(array_sphere, :p) |>x->reduce(hcat, x) |>x->collect(eachrow(x))
    scatter!(ax_points, x,y,z, color=D, markersize=30)
end



function update_pattern(res_pattern)
    # r = directivity_beta(res_pattern, θ_grid, ϕ_grid)
    r = 1 / 2(120pi) *(abs.(res_pattern.θ.coefs) .^2 .+  abs.(res_pattern.ϕ.coefs) .^2 )
    r = 10log10.(r) 
    x,y,z = begin
        local r_local =  r.+20
        replace(x->x < 0 ? 0 : x,  r)
        x =@. r_local*sin(θ_grid)cos(ϕ_grid)
        y =@. r_local*sin(θ_grid)sin(ϕ_grid)
        z =@. r_local*cos(θ_grid)
        [x,y,z]
    end
    [x,y,z,r]
end

lift(slider_θ.value, slider_ϕ.value) do θ,ϕ
    θ,ϕ = mod_angle_rad(deg2rad(θ), deg2rad(ϕ))
    @time res_pattern, Dir = cal_pattern(array_sphere, θ,ϕ, spin=true, optimal_directivity=:all);
    x_new,y_new,z_new,r_new = update_pattern(res_pattern)
    x[] = x_new    
    y[] = y_new    
    z[] = z_new    
    r[] = r_new    
    D[] = abs.(Dir).^2
    sleep(0)
end
fig



## ---------------Optimize---------------------
using Optim
