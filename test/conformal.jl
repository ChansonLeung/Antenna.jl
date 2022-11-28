using Polymake
using Antenna
using LinearAlgebra

point = Vector{Float64}

function get_center(points::Vector{point})
    sum(points)/length(points)
end

P = polytope.truncated_icosahedron()

vertices = P.VERTICES |> Matrix{Float64}|>eachrow.|> x->point(x[2:end])

# convert incidence Matrix to index list of points set
facets_index = P.VERTICES_IN_FACETS |> Matrix{Bool} |> eachrow .|>
  (row -> Iterators.filter(x->x[1],zip(row, 1:length(row)))) .|>
  row->map(i->i[2], row)

# get vector point
facets_point = map(x->vertices[x], facets_index)
# get face center point
face_center = map(get_center, facets_point)

#calculate the coordnate
coordinate = map(zip(face_center, 1:length(face_center))) do (point,index)
  # take point as z axis
  z = normalize(point)
  # take theta direction as x axis
  θ,ϕ,r = cart2sph(point...)
  x = [cos(θ)cos(ϕ), cos(θ)sin(ϕ), -sin(θ)]
  y = cross(z,x)
  ([x y z], index)
end

plot_coordinate(map(x->x[1], coordinate), face_center)