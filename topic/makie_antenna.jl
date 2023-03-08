using Antenna
using Makie

@recipe(Pattern3D) do scene
    Attributes(
        min_gain = -20
    )
end
@recipe(ArrayPoints) do scene
    Attributes(
        min_gain = -20
    )
end

@recipe(PatternUV) do scene
    Attributes(
        min_gain = -20
    )
end

function Makie.plot!(pattern3D_plot::Pattern3D{<:Tuple{anten_pattern}})
    pattern = pattern3D_plot[1]
    x = Observable(zeros(size(θ_grid)...))
    y = Observable(zeros(size(θ_grid)...))
    z = Observable(zeros(size(θ_grid)...))
    r = Observable(zeros(size(θ_grid)...))
    function update_plot(pattern)
        dir = 1 / 2(120pi) *(abs.(pattern.θ.coefs) .^2 .+  abs.(pattern.ϕ.coefs) .^2 )
        dir = 10log10.(dir) 
        local r_local =  dir.+ 20 |> r-> replace(x->x < 0 ? 0 : x,  r)
        x[] = @. r_local*sin(θ_grid)cos(ϕ_grid)
        y[] = @. r_local*sin(θ_grid)sin(ϕ_grid)
        z[] = @. r_local*cos(θ_grid)
        r[] = r_local
    end
    Makie.Observables.onany(update_plot, pattern)
    update_plot(pattern[])
    surface!(pattern3D_plot, x,y,z,color=r, colormap="lightrainbow")
end


function Makie.plot!(pattern3D_plot::Pattern3D{<:Tuple{anten_pattern}})
    pattern = pattern3D_plot[1]
    x = Observable(zeros(size(θ_grid)...))
    y = Observable(zeros(size(θ_grid)...))
    z = Observable(zeros(size(θ_grid)...))
    r = Observable(zeros(size(θ_grid)...))
    function update_plot(pattern)
        dir = 1 / 2(120pi) *(abs.(pattern.θ.coefs) .^2 .+  abs.(pattern.ϕ.coefs) .^2 )
        dir = 10log10.(dir) 
        local r_local =  dir.+ 20 |> r-> replace(x->x < 0 ? 0 : x,  r)
        x[] = @. r_local*sin(θ_grid)cos(ϕ_grid)
        y[] = @. r_local*sin(θ_grid)sin(ϕ_grid)
        z[] = @. r_local*cos(θ_grid)
        r[] = r_local
    end
    Makie.Observables.onany(update_plot, pattern)
    update_plot(pattern[])
    surface!(pattern3D_plot, x,y,z,color=r, colormap="lightrainbow")
end

function Makie.plot!(point_array::ArrayPoints{<:Tuple{Vector{anten_point}, Vector{Float64}}})
    array = point_array[1]
    input_D = point_array[2]
    x = Observable(zeros(size(array.val)...))
    y = Observable(zeros(size(array.val)...))
    z = Observable(zeros(size(array.val)...))
    c = Observable(zeros(size(input_D.val)...))

    function update_plot(array, input_D)
        x[] = getfield.(array, :p) .|> x->x[1]
        y[] = getfield.(array, :p) .|> x->x[2]
        z[] = getfield.(array, :p) .|> x->x[3]
        c[] = input_D
    end
    Makie.Observables.onany(update_plot, array, input_D)
    update_plot(array[], input_D[])
    scatter!(point_array, x,y,z,color=c)
end
