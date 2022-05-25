export
    scale2log,
    normlog,
    sph2cart,
    cart2sph,
    mod_angle,
    mod_angle_deg,
    mod_angle_rad


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
function mod_angle(θ, ϕ)
    sign = x -> x >= 0 ? 1 : -1
    item1 = -sign.(mod2pi.(θ) .- π)
    θ_ = mod1.(θ .* item1, π)
    θ_[θ.==0] .= 0
    ϕ_ = item1 .* ϕ .|> mod2pi
    (θ_, ϕ_)
end
function mod_angle(θ::Number, ϕ)
    sign = x -> x >= 0 ? 1 : -1
    item1 = -sign(mod2pi(θ) - π)
    θ_ = (θ==0 ?  0.0 : mod1(θ * item1, π))
    ϕ_ = (item1 .* ϕ) .|> mod2pi
    (θ_, ϕ_)
end

# convert to domain θ ∈[0°, 180°], ϕ∈ [-180°, 180°]
function mod_angle_deg(θ::Number, ϕ)
    θ_res = 0.
    ϕ_res = 0.
    invert_ϕ = false
    # in domain
    if θ >= 0 && θ <= 180
        θ_res = θ
    elseif θ < 0 && θ >= -180
        invert_ϕ = true 
        θ_res = -θ
    else
        throw("θ = $(θ) is out of range [-180,180]")
    end

    if ϕ >=-180 && ϕ <= 180
        ϕ_res = ϕ
    elseif ϕ > 180 && ϕ <=360
        ϕ_res = ϕ - 360
    else
        throw("ϕ = $(ϕ) is out of range [-180,360]")
    end

    if invert_ϕ
        if ϕ_res <= 0
            ϕ_res += 180
        else
            ϕ_res -= 180
        end
    end
    θ_res, ϕ_res
end

# not test
function mod_angle_rad(θ::Number, ϕ)
    θ_res = 0.
    ϕ_res = 0.
    invert_ϕ = false
    # in domain
    if θ >= 0 && θ <= pi
        θ_res = θ
    elseif θ < 0 && θ >= -pi
        invert_ϕ = true 
        θ_res = -θ
    else
        throw("θ = $(θ) is out of range [-pi, pi]")
    end

    if ϕ >=-pi && ϕ <= pi
        ϕ_res = ϕ
    elseif ϕ > pi && ϕ <=2pi
        ϕ_res = ϕ - 2pi
    else
        throw("ϕ = $(ϕ) is out of range [-pi, 2pi]")
    end

    if invert_ϕ
        if ϕ_res <= 0
            ϕ_res += pi
        else
            ϕ_res -= pi
        end
    end
    θ_res, ϕ_res
end