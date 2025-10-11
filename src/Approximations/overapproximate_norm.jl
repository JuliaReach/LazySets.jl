"""
    overapproximate_norm(Z::AbstractZonotope, [p]::Real=1)

Compute an upper bound on the ``p``-norm of a zonotopic set.

### Input

- `Z` -- Zonotopic set
- `p` -- (optional, default: `1`) norm

### Output

An upper bound on ``\\max_{x ∈ Z} ‖x‖_p``.

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

The problem ``\\max_{z ∈ Z} ‖z‖_1`` is NP-hard in general.

### Algorithm

This method computes a convex overapproximation of

```math
\(\max_{z∈Z}‖z‖₁\).
```
For each coordinate ``z_i ∈ [ℓ_i,u_i]``, the absolute value ``|z_i|`` is
replaced by its linear secant envelope:

```math
|z_i| ≤ a_i z_i + b_i, \quad
a_i = \frac{u_i + ℓ_i}{u_i - ℓ_i}, \quad
b_i = -\frac{2u_iℓ_i}{u_i - ℓ_i}.
```
Summing these bounds yields a linear function in ``z``, which is maximized
analytically over the zonotope as

```math
aᵀc + \sum_j |aᵀg_j| + \sum_i b_i.
```
This yields a linear-programming
relaxation of complexity \\(O(p·n)\\), where `p` is the number of generators and
`n` the ambient dimension.
"""
function _overapproximate_l1_norm(Z::AbstractZonotope{N}) where {N}
    lb = low(Z)
    ub = high(Z)
    n = dim(Z)

    S⁺ = Int[]
    S⁻ = Int[]
    S = Int[]
    @inbounds for i in 1:n
        if lb[i] ≥ 0
            push!(S⁺, i)
        elseif ub[i] ≤ 0
            push!(S⁻, i)
        else
            push!(S, i)
        end
    end

    w = ones(N, n)
    w[S⁻] .= N(-1)
    @. w[S] = (ub[S] + lb[S]) / (ub[S] - lb[S])

    const_term = -2 * sum(ub[S] .* lb[S] ./ (ub[S] .- lb[S]))

    # max_{x ∈ Z} wᵀx
    c = center(Z)
    G = genmat(Z)
    linear_max = dot(w, c) + sum(abs.(At_mul_B(G, w)))
    return linear_max + const_term
end

"""
    overapproximate_norm(MZ::MatrixZonotope, [p]::Real=Inf)

Compute an upper bound on the ``p``-norm of a matrix zonotope.

### Input

- `MZ` -- Matrix zonotope
- `p` -- (optional, default: `Inf`) norm

### Output

An upper bound on ``\\sup_{A ∈ \\mathcal{A} } ‖A‖_p ``.
"""
function overapproximate_norm(MZ::MatrixZonotope, p::Real=Inf)
    if p == 1
        return _rowwise_zonotope_norm(transpose(MZ), overapproximate_norm)
    elseif p == Inf
        return _rowwise_zonotope_norm(MZ, overapproximate_norm)
    else
        throw(ArgumentError("the norm for p=$p has not been implemented"))
    end
end

"""
    overapproximate_norm(MZP::MatrixZonotopeProduct, [p]::Real=Inf)

Compute an upper bound on the ``p``-norm of a product of matrix zonotopes.

### Input

- `MZP` -- Matrix zonotope product
- `p` -- (optional, default: `Inf`) norm

### Output

An upper bound on ``\\sup_{M ∈ \\mathcal{A*B*…} } \\|M\\|_p ``.
"""
function overapproximate_norm(MZP::MatrixZonotopeProduct, p::Real=Inf)
    return mapreduce(x -> overapproximate_norm(x, p), *, factors(MZP))
end
