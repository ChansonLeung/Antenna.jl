function Taylor_Chebyshev_I(n::Integer, R0::Number)::Function
    # Antenna Theory 7-26
    A_(R0) = 1 / π * acosh(R0)
    # Antenna Theory 7-28
    σ_(n, A) = n / √(A^2 + (n - 0.5)^2)
    # Antenna Theory 7-29
    μ_(n, σ, A) = π * σ * √(A^2 + (n - 0.5)^2)

    A = A_(R0)
    σ = σ_(n, A)
    μ(n) = μ_(n, σ, A)

    # Antenna Theory 7-30a
    SF(p) = begin
        item1 = factorial(n - 1)^2 / (factorial(n - 1 + p)factorial(n - 1 - p))
        item2 = prod([1 - (p * π / μ(i))^2 for i in 1:n-1])
        item1 * item2
    end
    # Antenna Theory 7-30
    I(l, z) = begin
        item1 = 1 / l
        item2 = sum(SF(p)cos(2π * p * z / l) for p in 1:n-1)
        item1 * (1 + 2 * item2)
    end
end