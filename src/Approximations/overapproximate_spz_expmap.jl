function _compute_series_coefficients(B, S0::S, k::Int) where {S}
    coeffs = Vector{S}(undef, k + 1)
    coeffs[1] = S0
    curr = S0
    # standard taylor: (1/i) * B * previous_term
    for i in 1:k
        curr = overapproximate(scale(1 / i, B) * curr, S)
        coeffs[i + 1] = curr
    end
    return coeffs
end

@inline function _horner_series(A, coeffs::Vector{S}) where {S}
    acc = coeffs[end]
    @inbounds for i in (length(coeffs) - 1):-1:1
        acc = exact_sum(overapproximate(A * acc, S), coeffs[i])
    end
    return acc
end

function _truncated_operator_series(MZP::MatrixZonotopeProduct,
                                    P::S, k::Int) where {S}
    facts = factors(MZP)
    B_last = facts[end]

    coeffs = _compute_series_coefficients(B_last, P, k)

    # fold with Horner's Method
    A_factors = reverse(view(facts, 1:(nfactors(MZP) - 1)))
    @inbounds for A in A_factors
        coeffs = [_horner_series(A, coeffs)]
    end
    return coeffs[1]
end

function _operator_series_remainder(Z::AbstractZonotope{N},
                                    _matnorm::Real, k::Int) where {N}
    ϵ = _matnorm / (k + 2)
    (ϵ >= 1) && @warn "Remainder series diverges; ε ≥ 1" ϵ k α

    n = dim(Z)
    E = IntervalMatrix(fill(IA.interval(N(-1.0), N(1.0)), n, n))
    factor = _matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

    R = overapproximate(scale(factor, E) * Z, Zonotope)
    return R
end

function _overapproximate_emz_generic(MZP::MAT, 
                                         P::S,
                                         k::Int,
                                         _matnorm::Real) where {S,MAT}
    require(@__MODULE__, :IntervalMatrices; fun_name="overapproximate")

    poly = _truncated_operator_series(MZP,P, k)

    Z = P isa AbstractZonotope ? P : overapproximate(P, Zonotope)
    rem = _operator_series_remainder(Z, _matnorm, k)

    return remove_redundant_generators(minkowski_sum(poly, rem))
end

"""
    overapproximate(em::ExponentialMap{N, S, MAT}, ::Type{<:S}, k::Int=2; matnorm=nothing) where {N, S, MAT}

Overapproximate the exponential map of a zonotopic set through a product of matrix zonotopes,
following Proposition 3 of [HuangLBS2025](@citet).

### Input

- `em`      -- exponential map of a zonotopic set (e.g., `Zonotope` or `SparsePolynomialZonotope`) 
               through a product of matrix zonotopes
- `Type`    -- target set type (matches the type of the set in `em`)
- `k`       -- (optional, default: `2`) order of the Taylor expansion
- `matnorm` -- (optional, default: `nothing`) pre-computed induced infinity-norm of the matrix zonotope product

### Output

A set of type `S` (the same type as the input set `set(em)`) representing the overapproximation.

### Notes

This function implements a Taylor-series-based overapproximation using Horner's method to minimize 
generator growth. It computes:

``\\bigoplus_{i=0}^k \\frac{A^i}{i!} S_i \\oplus \\mathcal{R}``

where ``S_i`` are coefficients determined by the expansion strategy and ``\\mathcal{R}`` is the 
Lagrange remainder bound.
"""
function overapproximate(em::ExponentialMap{N,S,MAT}, ::Type{<:S},
                         k::Int=2;
                         matnorm=nothing) where {N,S,MAT}
    _matnorm = isnothing(matnorm) ? N(overapproximate_norm(matrix(em).M, Inf)) : N(matnorm)

    # fallback to exact version if there are no generators (static matrix)
    MZP = matrix(em).M
    if MZP isa MatrixZonotope && isempty(generators(MZP))
        return exponential_map(center(MZP), set(em))
    end

    # homogeneous expansion (Prop 3 @HuangLBS2025)
    P = set(em)
    return _overapproximate_emz_generic(MZP, P, k, _matnorm)
end