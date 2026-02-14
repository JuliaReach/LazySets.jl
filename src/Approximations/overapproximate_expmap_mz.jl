function _compute_series_coefficients(B, S0::S,
                                      k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                              AbstractZonotope}}
    coeffs = Vector{S}(undef, k + 1)
    coeffs[1] = S0
    curr = S0

    @inbounds for i in 1:k
        curr = overapproximate(scale(1 / i, B) * curr, S)
        coeffs[i + 1] = curr
    end
    return coeffs
end

@inline function _horner_series(A::MatrixZonotope,
                                coeffs::Vector{S}) where {S<:Union{SparsePolynomialZonotope,
                                                                   AbstractZonotope}}
    acc = coeffs[end]
    @inbounds for i in (length(coeffs) - 1):-1:1
        acc = exact_sum(overapproximate(A * acc, S), coeffs[i])
    end
    return acc
end

function _truncated_operator_series(A::MatrixZonotope,
                                    P::S,
                                    k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                            AbstractZonotope}}
    return reduce(exact_sum, _compute_series_coefficients(A, P, k))
end

function _truncated_operator_series(MZP::MatrixZonotopeProduct,
                                    P::S,
                                    k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                            AbstractZonotope}}
    facts = factors(MZP)

    coeffs = _compute_series_coefficients(facts[end], P, k)

    # fold the previous factors with Horner's method
    for A in reverse(view(facts, 1:(nfactors(MZP) - 1)))
        coeffs = [_horner_series(A, coeffs)]
    end

    return coeffs[1]
end

function load_intervalmatrices_overapproximation_expmap_spz()
    return quote
        using .IntervalMatrices: IntervalMatrix, scale!

        function _operator_series_remainder(Z::AbstractZonotope{N},
                                            _matnorm::Real,
                                            k::Int) where {N}
            ϵ = _matnorm / (k + 2)
            if ϵ >= 1
                @warn "k should be chosen such that ϵ<1" ϵ
            end

            n = dim(Z)
            E = IntervalMatrix(fill(IA.interval(N(-1.0), N(1.0)), n, n))
            factor = _matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

            return overapproximate(scale!(E, factor) * Z, Zonotope)
        end
    end
end

function _overapproximate_emz_generic(MZP::MAT,
                                      P::S,
                                      k::Int,
                                      _matnorm::Real) where {S<:Union{SparsePolynomialZonotope,
                                                                      AbstractZonotope},
                                                             MAT<:AbstractMatrixZonotope}
    require(@__MODULE__, :IntervalMatrices; fun_name="overapproximate")

    tayexp = _truncated_operator_series(MZP, P, k)
    Z = P isa AbstractZonotope ? P : overapproximate(P, Zonotope)
    lagrem = _operator_series_remainder(Z, _matnorm, k)

    return remove_redundant_generators(minkowski_sum(tayexp, lagrem))
end

"""
    overapproximate(em::ExponentialMap{N,S,MAT},
                         ::Type{<:P},
                         k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N, S,
                                                                      P<:Union{SparsePolynomialZonotope,
                                                                               Zonotope},
                                                                      MAT<:AbstractMatrixZonotope{N}}

Overapproximate the exponential map of a set `X` of type `SparsePolynomialZonotope` or `Zonotope`
through a composition of matrix zonotopes, following Proposition 3 of [HuangLBS2025](@citet).

### Input

- `em`       -- an expontial map of set `X` through a product of matrix zonotopes
- `P`        -- target type
- `k`        -- (default: `2`) the order of the taylor expansion
- `matnorm`  -- (Optional, default: `nothing`) Pre-computed induced ``\\infty``-norm of the matrix zonotope

### Output

A set overapproximating the exponential map.
"""
function overapproximate(em::ExponentialMap{N,S,MAT},
                         ::Type{<:P},
                         k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N, S,
                                                                      P<:Union{SparsePolynomialZonotope,
                                                                               Zonotope},
                                                                      MAT<:AbstractMatrixZonotope{N}}
    _matnorm = isnothing(matnorm) ? N(overapproximate_norm(matrix(em).M, Inf)) : N(matnorm)

    MZP = matrix(em).M
    if MZP isa MatrixZonotope && isempty(generators(MZP))
        return exponential_map(center(MZP), set(em))
    end

    return _overapproximate_emz_generic(MZP, set(em), k, _matnorm)
end

function overapproximate(em::ExponentialMap{N,S,MAT},
                         ::Type{<:P},
                         k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N, S,
                                                                      P<:Union{SparsePolynomialZonotope,
                                                                               Zonotope},
                                                                      MAT<:SparseMatrixExp{N}}
    return linear_map(matrix(em), set(em))
end
