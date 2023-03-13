using Antenna
using Makie
using Meshes
    


Makie.@recipe(Pattern3D) do scene
    Makie.Attributes(
        min_gain = -20
    )
end
Makie.@recipe(ArrayPoints) do scene
    Makie.Attributes(
        min_gain = -20
    )
end

Makie.@recipe(PatternUV) do scene
    Makie.Attributes(
        min_gain = -20
    )
end

function Makie.plot!(pattern3D_plot::Pattern3D{<:Tuple{Matrix{Float64}}})
    pattern = pattern3D_plot[1]
    x = Observable(zeros(size(θ_grid)...))
    y = Observable(zeros(size(θ_grid)...))
    z = Observable(zeros(size(θ_grid)...))
    r = Observable(zeros(size(θ_grid)...))
    function update_plot(pattern)
        dir = pattern
        dir = 10log10.(dir) 
        local r_local =  dir.+ 20 |> r-> replace(x->x < 0 ? 0 : x,  r)
        x[] = @. r_local*sin(θ_grid)cos(ϕ_grid)
        y[] = @. r_local*sin(θ_grid)sin(ϕ_grid)
        z[] = @. r_local*cos(θ_grid)
        r[] = r_local
    end
    Makie.Observables.onany(update_plot, pattern)
    update_plot(pattern[])
    Makie.surface!(pattern3D_plot, x,y,z,color=r, colormap="lightrainbow")
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
    Makie.surface!(pattern3D_plot, x,y,z,color=r, colormap="lightrainbow")
end


# function Makie.plot!(point_array::ArrayPoints{<:Tuple{Vector{anten_point}, Vector{Float64}}})
#     array = point_array[1]
#     input_D = point_array[2]
#     x = Observable(zeros(size(array.val)...))
#     y = Observable(zeros(size(array.val)...))
#     z = Observable(zeros(size(array.val)...))
#     c = Observable(zeros(size(input_D.val)...))

#     function update_plot(array, input_D)
#         x[] = getfield.(array, :p) .|> x->x[1]
#         y[] = getfield.(array, :p) .|> x->x[2]
#         z[] = getfield.(array, :p) .|> x->x[3]
#         c[] = input_D
#     end
#     Makie.Observables.onany(update_plot, array, input_D)
#     update_plot(array[], input_D[])
#     Makie.scatter!(point_array, x,y,z,color=c)
# end

# function Makie.plot!(point_array::ArrayPoints{<:Tuple{Vector{Vector{Float64}}, Vector{Float64}}})
#     array = point_array[1]
#     color = point_array[2]
#     x = Observable(zeros(size(array.val)...))
#     y = Observable(zeros(size(array.val)...))
#     z = Observable(zeros(size(array.val)...))
#     c = Observable(zeros(size(color.val)...))

#     function update_plot(array,color)
#         x[] = array .|> x->x[1]
#         y[] = array .|> x->x[2]
#         z[] = array .|> x->x[3]
#         c[] = color
#     end
#     Makie.Observables.onany(update_plot, array, color)
#     update_plot(array[], color[])
#     Makie.scatter!(point_array, x,y,z, color=c)
# end


Makie.convert_arguments(plottype::Type{<:ArrayPoints}, vec_p::Vector{Vector{Float64}}) = 
    convert_arguments(plottype, vec_p, zeros(size(vec_p)))

Makie.convert_arguments(plottype::Type{<:ArrayPoints}, vec_p::Vector{Meshes.Point3}) = 
    convert_arguments(plottype, map(x->[x.coords...], vec_p))

## dispatch scatter
Makie.convert_arguments(plottype::Type{<:Scatter}, vec_p::Vector{Vector{Float64}}) = 
    convert_arguments(plottype, map(x->Makie.Point3(x), vec_p))
Makie.convert_arguments(plottype::Type{<:Scatter}, vec_p::Vector{anten_point}) = 
    convert_arguments(plottype, map(x->Makie.Point3(x.p), vec_p))
Makie.convert_arguments(plottype::Type{<:Scatter}, vec_p::Vector{Meshes.Point3}) = 
    convert_arguments(plottype, map(x->Makie.Point3(x.coords...), vec_p))
# Makie.convert_arguments(plottype::Type{<:ArrayPoints}, vec_p::Vector{anten_point}) = 
#     convert_arguments(plottype, map(x->x.p, vec_p), map(x->x.coefs, vec_p))