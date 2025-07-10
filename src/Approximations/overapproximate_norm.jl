"""
    overapproximate_norm(Z::AbstractZonotope, [p]::Real=1)

Compute an upper bound on the ``p``-norm of a zonotopic set.

### Input

- `Z` -- Zonotopic set
- `p` -- (optional, default: `1`) norm

### Output

An upper bound on ``\\max_{x ∈ Z} \\|x\\|_p``.

### Notes

The norm of a set is the norm of the enclosing ball (of the given ``p``-norm) of
minimal volume that is centered in the origin.
"""
function overapproximate_norm(Z::AbstractZonotope, p::Real=1)
    if p == 1
        return _overapproximate_l1_norm(Z)
    else
        return norm(Z, p)
    end
end

"""
    _overapproximate_l1_norm(Z::AbstractZonotope{N}) where {N}

Compute an upper bound on the 1-norm of a zonotopic set.

### Notes

The problem ``\\max_{z ∈ Z} \\|z\\|_1`` is NP-hard in general.

### Algorithm

This method computes a convex relaxation by combining the MILP formulation
described in [JordanD21; Theorem 5](@citet) with the convex-hull construction in
the supplement section §C.3. We replace each coordinate's absolute value

```math
|z_i|\\quad z_i ∈ [ℓ_i, u_i]
```

by its convex-hull secant upper envelope: the line through the two points
``((ℓ_i,|ℓ_i|)`` and ``((u_i,|u_i|)``. This yields a linear-programming
relaxation of complexity \\(O(p·n)\\), where `p` is the number of generators and
`n` the ambient dimension.
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
