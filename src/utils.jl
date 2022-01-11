module utils

export 
    scale2log,
    normlog,
    sph2cart,
    cart2sph

# operation
scale2log(x) = 20log10(x)
normlog(vec) = vec .- maximum(vec)
sph2cart(θ, ϕ, r) = begin
    (
        x = r * sin(θ)cos(ϕ),
        y = r * sin(θ)sin(ϕ),
        z = r * cos(θ)
    )
end
cart2sph(x, y, z) = begin
    # XXX r are not calcuated, may cause problem
    r = sqrt(x^2 + y^2 + z^2)
    θ = acos(z / r)
    ϕ = atan(y, x)
    (
        θ = θ,
        ϕ = ϕ,
        r = r
    )
end

end