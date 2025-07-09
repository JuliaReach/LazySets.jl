"""
    overapproximate_norm(Z::AbstractZonotope, p::Real=1)

Compute an upper bound on the ``\\ell_p`` norm over a zonotope `Z`.

### Input

- `Z` -- A zonotopic set
- `p` -- The ``ℓ_p`` norm to overapproximate

### Output

- An upper bound on ``\\max_{z ∈ Z} \\|z\\|_p``.
"""
function overapproximate_norm(Z::AbstractZonotope, p::Real=1)
    if p == 1
        return _overapproximate_l1_norm(Z)
    else
        return norm(Z, p)
    end
end

"""
    _overapproximate_l1_norm(Z::AbstractZonotope)

Compute an upper bound on the ``\\ell_1`` norm over a zonotope `Z`.

### Notes

The problem ``\\max_{z ∈ Z} \\|z\\|_1`` is NP-hard in general. This function computes a convex relaxation 
by combining the MILP formulation described in [Jordan2021; Theorem 5](@citet) with the convex-hull construction
in their supplement section §C.3.

### Algorithm 

We replace each coordinate's absolute value
``
|z_i|\\quad z_i ∈ [ℓ_i, u_i]
``
by its convex-hull secant upper envelope: the line through
the two points ``((ℓ_i,|ℓ_i|)`` and ``((u_i,|u_i|)``.  This yields a
linear-programming relaxation of complexity \\(O(n·d)\\), where `n` is the
number of generators and `d` the ambient dimension.
"""
function _overapproximate_l1_norm(Z::AbstractZonotope{N}) where {N}
    lb = low(Z)
    ub = high(Z)

    c = center(Z)
    G = genmat(Z)

    S⁺ = lb .> 0
    S⁻ = ub .< 0
    S = .!(S⁺ .| S⁻)

    w = zeros(N, length(c))
    w[S⁺] .= N(1)
    w[S⁻] .= N(-1)
    @. w[S] = (ub[S] + lb[S]) / (ub[S] - lb[S])

    const_term = -2 * sum(ub[S] .* lb[S] ./ (ub[S] .- lb[S]))
    # max_{x ∈ Z} wᵀx
    linear_max = dot(w, c) + sum(abs.(At_mul_B(G, w)))
    return linear_max + const_term
end
