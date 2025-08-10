"""
    _compute_inner_powers(B::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}

Compute the first `k` powers of the matrix zonotope `B` applied to the set `P`.

This function returns a vector of overapproximated sets of the form:

```math
    (1/i!) * \\mathcal{B}^i * P
```
for `i = 0` to `k`.
"""
function _compute_inner_powers(B::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    powers = Vector{S}(undef, k + 1)
    powers[1] = P

    @inbounds for i in 1:k
        term = overapproximate(scale(1 / i, B) * powers[i], S)
        powers[i + 1] = term
    end
    return powers
end

"""
    _compute_outer_powers(A::MatrixZonotope, in_powers::Vector{S},
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}

Apply `A` repeatedly to each element in `in_powers`, approximating:

```math
A^i * (B^i * P)
```

for `i = 0` to `k`.
"""
function _compute_outer_powers(A::MatrixZonotope, in_powers::Vector{S},
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    out_powers = similar(in_powers)
    out_powers[1] = in_powers[1]

    @inbounds for i in 2:(k + 1)
        term = in_powers[i]
        for _ in 1:(i - 1)
            term = overapproximate(A * term, S)
        end
        out_powers[i] = term
    end
    return out_powers
end

"""
    taylor_expmap_truncation(A::MatrixZonotope, P::S, k::Int)
        where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}

Compute the k-th order truncated Taylor expansion of the exponential map of a matrix zonotope

### Input

- `A` -- a matrix zonotope
- `P` -- a (potentially polynomial) zonotopic set
- `k` -- the order of the Taylor expansion

### Output

A (polynomial) zonotopic set representing the k-th order truncated Taylor expansion.

### Algorithm

This function computes the approximation:

```math
\\displaystyle\\boxplus_{i=0}^k \\frac{\\mathcal{A}^i }{i!} X
```
"""
function taylor_expmap_truncation(A::MatrixZonotope, P::S,
                                  k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                          AbstractZonotope}}
    inner = _compute_inner_powers(A, P, k)
    res = reduce(exact_sum, inner)
    return res
end

"""
    taylor_expmap_truncation(MZP::MatrixZonotopeProduct, P::S, k::Int)
        where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}

Compute the k-th order truncated Taylor expansion of the exponential map of a matrix zonotope product

### Input

- `MZP` -- a matrix zonotope product
- `P` -- a (potentially polynomial) zonotopic set
- `k` -- the order of the Taylor expansion

### Output

A (potentially polynomial) zonotopic set representing the k-th order truncated Taylor expansion.

### Algorithm

This function computes the approximation:

```math
\\displaystyle\\boxplus_{i=0}^k \\frac{\\mathcal{A}^i \\mathcal{B}^i}{i!} P
```
"""
function taylor_expmap_truncation(MZP::MatrixZonotopeProduct, P::S,
                                  k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                          AbstractZonotope}}
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
            taylor_expmap_remainder(Z::AbstractZonotope{N}, matnorm::Real, Int) where {N}

        Overapproximate the Lagrange remainder term of the k-th order truncated Taylor expansion
        of the exponential map of a matrix zonotope applied to a zonotopic set.

        ### Input

        - `P` -- a zonotopic set
        - `matnorm` -- an upper bound on the norm of the matrix zonotope
        - `k` -- the order of the Taylor expansion

        ### Output

        A zonotope over-approximating the remainder term of the Taylor expansion.
        """
        function taylor_expmap_remainder(Z::AbstractZonotope{N}, matnorm::Real, k::Int) where {N}
            n = dim(Z)
            ϵ = matnorm / (k + 2)

            E = IntervalMatrix(fill(IA.interval(N(-1.0), N(1.0)), n, n))
            E *= matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

            res = overapproximate(E * Z, Zonotope)
            return res
        end

        """
            overapproximate(expA::MatrixZonotopeExp{N,T}, ::Type{<:MatrixZonotope},
                                 k::Int=2) where {N,T<:AbstractMatrixZonotope{N}}
        
        Overapproximate the matrix zonotope exponential ``exp(\\mathacal{A})``

        ### Input

        - `expA` -- A `MatrixZonotopeExp`
        - `MatrixZonotope` -- target type
        - `k` -- (default: `2`) the order of the taylor expansion

        ### Output 
        
        A matrix zonotope overapproximating the matrix zonotope exponential

        ### Algorithm
        
        The algorithm follows [AlthoffKS11](@citet). 
        The expansion

        ```math 
        exp(\\mathcal{A}) ⊆ \\sum_i^k \\frac{\\mathcal{A}^i}{i!} + E_k
        ```

        is computed by overappraximating the matrix zonotope powers ``A^i`` 
        for ``i=0, \\dots, k``. 
        The remainder term ``E_k`` is computed through interval arithmetic 
        following Proposition 4.1 
        """
        function overapproximate(expA::MatrixZonotopeExp{N,T}, ::Type{<:MatrixZonotope},
                                 k::Int=2) where {N,T<:AbstractMatrixZonotope{N}}
            # overapproximate the exponent A*B*...*D
            MZP = MatrixZonotopeProduct(expA.M)
            X = overapproximate(MZP, MatrixZonotope)

            # compute the taylor expansion 
            powers = Vector{typeof(X)}(undef, k)
            powers[1] = X
            @inbounds for i in 2:k
                term = overapproximate(X * powers[i - 1], MatrixZonotope)
                powers[i] = scale(1 / i, term)
            end
            W = reduce(minkowski_sum, powers)
            W = MatrixZonotope(center(W) + Matrix{N}(I, size(W)), generators(W))

            # overapproximate mat zon by interval matrix and overapproximate remainder
            A = overapproximate(X, IntervalMatrix)
            E = IntervalMatrices._exp_remainder(A, N(1), k)
            
            res = minkowski_sum(W, convert(MatrixZonotope, E))
            #TODO change to remove_redundant_generators(res) after closing #3999 
            return res
        end
    end
end

function _overapproximate_emz(em::ExponentialMap{N,S,MAT}, k::Int,
                              matnorm::Real) where {N,
                                                    S<:Union{SparsePolynomialZonotope,
                                                             AbstractZonotope},
                                                    MAT<:AbstractMatrixZonotope{N}}
    require(@__MODULE__, :IntervalMatrices; fun_name="overapproximate")

    MZP = matrix(em).M
    P = set(em)

    if MZP isa MatrixZonotope && isempty(generators(MZP))
        return exponential_map(center(MZP), P)
    end

    ϵ = matnorm / (k + 2)
    if ϵ >= 1
        @warn "k should be chosen such that ϵ<1 " ϵ
    end

    tayexp = taylor_expmap_truncation(MZP, P, k)
    Z = P isa AbstractZonotope ? P : overapproximate(P, Zonotope)
    lagrem = taylor_expmap_remainder(Z, matnorm, k)
    P_approx = minkowski_sum(tayexp, lagrem)

    return remove_redundant_generators(P_approx)
end

"""
    overapproximate(em::ExponentialMap{N,S,MAT},
                         ::Type{<:SparsePolynomialZonotope},
                         k::Int=2;
                         [matnorm]::Union{Real,Nothing}=nothing) where {N,
                                                                      S<:SparsePolynomialZonotope,
                                                                      MAT<:AbstractMatrixZonotope{N}}

Overapproximate the exponential map of a sparse polynomial zonotope through a composition of matrix
zonotopes, following Proposition 3 of [HuangLBS2025](@citet).

### Input

- `em`      -- an expontial map of a sparse polynomial zonotope through a product of matrix zonotopes
- `SparsePolynomialZonotope` -- target type
- `k`       -- (default: `2`) the order of the taylor expansion
- `matnorm` -- (Optional, default: `nothing`) Pre-computed induced ``\\infty``-norm of the matrix zonotope

### Output

A sparse polynomial zonotope overapproximating the exponential map.
"""
function overapproximate(em::ExponentialMap{N,S,MAT},
                         ::Type{<:SparsePolynomialZonotope},
                         k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N,
                                                                      S<:SparsePolynomialZonotope,
                                                                      MAT<:AbstractMatrixZonotope{N}}
    _matnorm = isnothing(matnorm) ? N(overapproximate_norm(matrix(em).M, Inf)) : N(matnorm)
    return _overapproximate_emz(em, k, _matnorm)
end

"""
    overapproximate(em::ExponentialMap{N,S,MAT}, ::Type{<:Zonotope},
                         k::Int=2;
                         [matnorm]::Union{Real,Nothing}=nothing) where {N,S<:AbstractZonotope,
                                                                      MAT<:AbstractMatrixZonotope{N}}

Overapproximate the exponential map of a zonotope through a composition of matrix
zonotopes, following Proposition 3 of [HuangLBS2025](@citet).

### Input

- `em`       -- an expontial map of a zonotope through a product of matrix zonotopes
- `Zonotope` -- target type
- `k`        -- (default: `2`) the order of the taylor expansion
- `matnorm`  -- (Optional, default: `nothing`) Pre-computed induced ``\\infty``-norm of the matrix zonotope

### Output

A zonotope overapproximating the exponential map.
"""
function overapproximate(em::ExponentialMap{N,S,MAT}, ::Type{<:Zonotope},
                         k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N,S<:AbstractZonotope,
                                                                      MAT<:AbstractMatrixZonotope{N}}
    _matnorm = isnothing(matnorm) ? N(overapproximate_norm(matrix(em).M, Inf)) : N(matnorm)
    return _overapproximate_emz(em, k, _matnorm)
end

function overapproximate(em::ExponentialMap{N,S,MAT}, ::Type{<:SparsePolynomialZonotope}, k::Int=2;
                         matnorm::Union{Real,Nothing}=nothing) where {N,S<:SparsePolynomialZonotope,
                                                                      MAT<:SparseMatrixExp{N}}
    return linear_map(matrix(em), set(em))
end
