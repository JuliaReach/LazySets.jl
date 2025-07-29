"""
    _compute_inner_powers(A::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    
Compute the first `k` inner powers of the matrix zonotope `A` applied to the set `P`.

This function returns a vector of overapproximated sets of the form:

```math
    (1/n!) * \\mathcal{A}^n * P
```
for `n = 0` to `k`.
"""
function _compute_inner_powers(A::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    T = typeof(P)
    powers = Vector{T}(undef, k + 1)
    invfact = 1
    powers[1] = P

    @inbounds for n in 1:k
        invfact /= n
        term = overapproximate(scale(invfact, A) * powers[n], T)
        powers[n + 1] = term
    end
    return powers
end

"""
    _compute_outer_powers(B::MatrixZonotope, in_powers::Vector{S},
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}

Apply `B` repeatedly to each element in `in_powers`, approximating:

```math
B^n * (A^n * P)
```

for `n = 0` to `k`.
"""
function _compute_outer_powers(B::MatrixZonotope, in_powers::Vector{S},
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    out_powers = similar(in_powers)
    out_powers[1] = in_powers[1]

    @inbounds for n in 2:(k + 1)
        term = in_powers[n]
        for _ in 1:n
            term = overapproximate(B * term, S)
        end
        out_powers[n] = term
    end
    return out_powers
end

"""
    taylor_expmap_truncation(A::MatrixZonotope, P::S, k::Int) where S

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
function taylor_expmap_truncation(A::MatrixZonotope, P::S,
                        k::Int) where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}
    if isempty(generators(A))
        res = linear_map(exp(center(A)), P)
    else
        inner = _compute_inner_powers(A, P, k)
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
    for A in reverse(view(factors(MZP), 1:(nfactors(MZP) - 1)))
        terms = _compute_outer_powers(A, terms, k)
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

            tayexp = taylor_expmap_truncation(MZP, P, k)
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
