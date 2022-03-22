export
    scale2log,
    normlog,
    sph2cart,
    cart2sph,
    mod_angle

# operation
scale2log(x) = 20log10(x)
normlog(vec) = vec .- maximum(vec)
sph2cart(θ, ϕ, r) = begin
    (
        x=r * sin(θ)cos(ϕ),
        y=r * sin(θ)sin(ϕ),
        z=r * cos(θ)
    )
end
cart2sph(x, y, z) = begin
    # XXX r are not calcuated, may cause problem
    r = sqrt(x^2 + y^2 + z^2)
    θ = acos(z / r)
    ϕ = atan(y, x)
    (
        θ=θ,
        ϕ=ϕ,
        r=r
    )
end
# convert to domain θ ∈[0°, 180°], ϕ∈ [-180°, 180°]
# bad implementation
mod_angle = (θ, ϕ) -> begin
    sign = x -> x >= 0 ? 1 : -1
    item1 = -sign.(mod2pi.(θ) .- π)
    θ_ = mod1.(θ .* item1, π)
    θ_[θ.==0] .= 0
    ϕ_ = item1 .* ϕ .|> mod2pi
    (θ_, ϕ_)
end