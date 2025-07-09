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
using a formulation described in [Jordan2021; Theorem 5](@citet), which fits a tight upper envelope to the 
absolute value operator. For a zonotope ``Z ∈ ℝ^d`` with `n` generators, this function
has time complexity ``\\mathcal{O}(n ⋅ d)`` compared to ``\\mathcal{O}(2ⁿ ⋅ d)`` for the exact function.
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
    w[S⁺] .= 1.0
    w[S⁻] .= -1.0
    @. w[S] = (ub[S] + lb[S]) / (ub[S] - lb[S])

    const_term = -2 * sum(ub[S] .* lb[S] ./ (ub[S] .- lb[S]))
    # max_{x ∈ Z} wᵀx
    linear_max = dot(w, c) + sum(abs.(At_mul_B(G, w)))
    return linear_max + const_term
end
