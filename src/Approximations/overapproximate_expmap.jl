function _compute_inner_powers(MZ::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    T = typeof(P)
    powers = Vector{T}(undef, k + 1)
    invfact = 1
    powers[1] = P
    
    @inbounds for i in 1:k
        invfact /= i
        term = overapproximate(scale(invfact, MZ) * powers[i], T)
        powers[i + 1] = term
    end
    return powers
end

function _compute_outer_powers(MZ::MatrixZonotope, in_powers::Vector{S},
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    out_powers = similar(in_powers)
    out_powers[1] = in_powers[1]

    @inbounds for i in 2:(k + 1)
        term = in_powers[i]
        for _ in 1:i
            term = overapproximate(MZ * term, S)
        end
        out_powers[i] = term
    end
    return out_powers
end

"""
    taylor_expmap_truncation(MZ::MatrixZonotope, P::S, k::Int) where S

Compute the k-th order truncated Taylor expansion of the exponential map of a matrix zonotope

### Input

- `A` -- a matrix zonotope
- `P` -- a zonotopic set
- `k` -- the order of the Taylor expansion

### Output

A zonotopic set representing the k-th order truncated Taylor expansion.

### Algorithm

This function computes the approximation:

```math
\\displaystyle\\boxplus_{i=0}^k \\frac{\\mathcal{A}^i }{i!} X
```
"""
function taylor_expmap_truncation(MZ::MatrixZonotope, P::S,
                        k::Int) where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}
    if isempty(generators(MZ))
        res = linear_map(exp(center(MZ)), P)
    else
        inner = _compute_inner_powers(MZ, P, k)
        res = reduce(exact_sum, inner)
    end

    return res
end

"""
    taylor_expmap_truncation(MZP::MatrixZonotopeProduct, P::S, k::Int) where S

Compute the k-th order truncated Taylor expansion of the exponential map of a matrix zonotope product

### Input

- `MZP` -- a matrix zonotope product
- `P` -- a zonotopic set
- `k` -- the order of the Taylor expansion

### Output

A zonotopic set representing the k-th order truncated Taylor expansion.

### Algorithm

This function computes the approximation:

```math
\\displaystyle\\boxplus_{i=0}^k \\frac{\\mathcal{A}^i \\mathcal{B}^i}{i!} X
```
"""
function taylor_expmap_truncation(MZP::MatrixZonotopeProduct, P::S,
                        k::Int) where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}
    # inner powers on the last factor
    last_factor = factors(MZP)[end]
    terms = _compute_inner_powers(last_factor, P, k)

    # propagate outwards through the preceding factors
    for MZ in reverse(view(factors(MZP), 1:(nfactors(MZP) - 1)))
        terms = _compute_outer_powers(MZ, terms, k)
    end

    res = reduce(exact_sum, terms)
    return res
end

function load_intervalmatrices_overapproximation_expmap()
    return quote
        using .IntervalMatrices: IntervalMatrix

        """
            taylor_expmap_remainder(P::S, matnorm::Real, k::Int) where S

        Compute the Lagrange remainder term of the k-th order truncated Taylor expansion 
        of the exponential map of a matrix zonotope applied to a zonotopic set.

        ### Input

        - `P` -- a zonotopic set
        - `matnorm` -- an upper bound on the norm of the matrix zonotope
        - `k` -- the order of the Taylor expansion

        ### Output

        A zonotope over-approximating the remainder term of the Taylor expansion.
        """
        function taylor_expmap_remainder(P::S, matnorm::Real, k::Int) where{S<:Union{SparsePolynomialZonotope,
                                                                AbstractZonotope}}
            n = dim(P)
            N = eltype(P)

            ϵ = matnorm / (k + 2)
            ε = IntervalMatrix(fill(IA.interval(N(-1.0), N(1.0)), n, n))
            ε *= matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

            Zp = overapproximate(P, Zonotope)
            rhs = overapproximate(ε * Zp, Zonotope)
            return rhs
        end

        function _overapproximate_emz(em::ExponentialMap{N,S,MAT},
                                      k::Int=2) where {N,S<:Union{SparsePolynomialZonotope,
                                                                AbstractZonotope},
                                                       MAT<:AbstractMatrixZonotope{N}}
            MZP = matrix(em).M
            P = set(em)

            matnorm = N(overapproximate_norm(MZP, Inf))
            ϵ = matnorm / (k + 2)
            if ϵ > 1
                @warn "k should be chosen such that ϵ<1 " ϵ
            end

            tayexp = _taylor_expmap(MZP, P, k)
            lagrem = lagrange_remainder(P, matnorm, k)
            P_approx = minkowski_sum(tayexp, lagrem)

            return remove_redundant_generators(P_approx)
        end

        """
        	overapproximate(em::ExponentialMap{N,S,MAT},
                                 T::Type{<:SparsePolynomialZonotope},
                                 k::Int=2) where {N,S<:SparsePolynomialZonotope,
                                                  MAT<:AbstractMatrixZonotope{N}}

        Overapproximate the exponential map of a sparse polynomial zonotope through a composition of matrix 
        zonotopes, following Proposition 1 of [HuangLBS2025](@citet).

        ### Input

        - `em` -- an expontial map of a sparse polynomial zonotope through a product of matrix zonotopes

        ### Output

        A sparse polynomial zonotope overapproximating the exponential map.
        """
        function overapproximate(em::ExponentialMap{N,S,MAT},
                                 T::Type{<:SparsePolynomialZonotope},
                                 k::Int=2) where {N,S<:SparsePolynomialZonotope,
                                                  MAT<:AbstractMatrixZonotope{N}}
            return _overapproximate_emz(em, k)
        end

        """
        	overapproximate(em::ExponentialMap{N,S,MAT}, T::Type{<:Zonotope},
                                 k::Int=2) where {N,S<:AbstractZonotope,
                                                  MAT<:AbstractMatrixZonotope{N}}

        Overapproximate the exponential map of a zonotope through a composition of matrix 
        zonotopes, following Proposition 1 of [HuangLBS2025](@citet).

        ### Input

        - `em` -- an expontial map of a zonotope through a product of matrix zonotopes

        ### Output

        A zonotope overapproximating the exponential map.
        """
        function overapproximate(em::ExponentialMap{N,S,MAT}, T::Type{<:Zonotope},
                                 k::Int=2) where {N,S<:AbstractZonotope,
                                                  MAT<:AbstractMatrixZonotope{N}}
            return _overapproximate_emz(em, k)
        end
    end
end
