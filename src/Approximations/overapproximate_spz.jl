"""
	overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:SparsePolynomialZonotope}) where {N,S<:SparsePolynomialZonotope{N},
                                                                  NM,
                                                                  MAT<:MatrixZonotope{NM}}

Overapproximate the linear map of a sparse polynomial zonotope through a matrix zonotope,
following Proposition 2 of [HuangLBS2025](@citet).

### Input

- `lm` -- a linear map of a sparse polynomial zonotope through a matrix zonotope

### Output

A sparse polynomial zonotope overapproximating the linear map.

"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:SparsePolynomialZonotope}) where {N,S<:SparsePolynomialZonotope{N},
                                                                  NM,
                                                                  MAT<:MatrixZonotope{NM}}
    MZ = matrix(lm)
    P = set(lm)
    T = promote_type(N, NM)

    m, n = size(MZ)
    w = ngens(MZ)
    h = ngens_dep(P)
    q = ngens_indep(P)

    c = center(MZ) * center(P)

    # matrix of independent generators
    Gi = Matrix{T}(undef, m, q * (w + 1))
    Gi[:, 1:q] = center(MZ) * genmat_indep(P)

    # compute matrix of dependent generators
    G = Matrix{T}(undef, m, h + w + h * w)
    G[:, 1:h] = center(MZ) * genmat_dep(P)

    # loop to populate G and Gi
    @inbounds for (i, A) in enumerate(generators(MZ))
        G[:, h + i] = A * center(P)
        G[:, (h + w + (i - 1) * h + 1):(h + w + i * h)] = A * genmat_dep(P)
        Gi[:, (q * i + 1):(q * (i + 1))] = A * genmat_indep(P)
    end

    # compute exponent
    Imat = Matrix{Int}(I, w, w)
    Ê₁, Ê₂, idx = merge_id(indexvector(P), indexvector(MZ), expmat(P), Imat)
    pₖ = size(Ê₁, 1)
    E = Matrix{Int}(undef, pₖ, h + w + h * w)
    E[:, 1:h] = Ê₁
    E[:, (h + 1):(h + w)] = Ê₂
    ones_h = ones(Int, 1, h)
    @inbounds for l in 1:w
        cstart = (h + w) + (l - 1) * h + 1
        cend = (h + w) + l * h
        col = @view Ê₂[:, l]
        E[:, cstart:cend] = col * ones_h .+ Ê₁
    end    

    return SparsePolynomialZonotope(c, G, Gi, E, idx)
end

"""
	function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                  MAT<:MatrixZonotope{NM}}

Overapproximate the linear map of a zonotope through a matrix zonotope,
following a modification of Proposition 1 of [HuangLBS2025](@citet).

### Input

- `lm` -- a linear map of a zonotope through a matrix zonotope

### Output

A zonotope overapproximating the linear map.

"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                  MAT<:MatrixZonotope{NM}}
    MZ = matrix(lm)
    Z = set(lm)
    T = promote_type(N, NM)

    m, n = size(MZ)
    w = ngens(MZ)
    h = ngens(Z)

    c = center(MZ) * center(Z)

    # compute matrix of dependent generators
    G = Matrix{T}(undef, m, h + w + h * w)
    G[:, 1:h] = center(MZ) * genmat(Z)
    @inbounds for (i, A) in enumerate(generators(MZ))
        G[:, h + i] = A * center(Z)
        G[:, (h + w + (i - 1) * h + 1):(h + w + i * h)] = A * genmat(Z)
    end

    return Zonotope(c, G)
end

function _overapproximate_lmzp(lm::LinearMap{N,S,NM,MAT}) where {N,S<:AbstractZonotope{N},NM,
                                                                 MAT<:MatrixZonotopeProduct{NM}}
    MZs = factors(matrix(lm))
    P = set(lm)

    # apply overapproximation from innermost to outermost
    reduced = foldr((A, acc) -> overapproximate(A * acc, U), MZP_factors; init=P)
    return reduced
end

"""
    overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                    MAT<:MatrixZonotopeProduct{NM}}

Overapproximate the linear map of a zonotope through a product of matrix zonotopes,
by recursively applying the overapproximation rule from the inside out.

### Input

- `lm` -- a linear map of a zonotope through a `MatrixZonotopeProduct`
- `U` -- the target overapproximation type

### Output

An overapproximation of the linear map as a zonotope.
"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                    MAT<:MatrixZonotopeProduct{NM}}
    return _overapproximate_lmzp(lm)
end

"""
    overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:SparsePolynomialZonotope}) where {N,
                                                                    S<:SparsePolynomialZonotope{N},
                                                                    NM,
                                                                    MAT<:MatrixZonotopeProduct{NM}}

Overapproximate the linear map of a sparse polynomial zonotope through a product of matrix zonotopes,
by recursively applying the overapproximation rule from the inside out.

### Input

- `lm` -- a linear map of a sparse polynomial zonotope through a `MatrixZonotopeProduct`
- `U` -- the target overapproximation type

### Output

An overapproximation of the linear map as a sparse polynomial zonotope,
"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:SparsePolynomialZonotope}) where {N,
                                                                    S<:SparsePolynomialZonotope{N},
                                                                    NM,
                                                                    MAT<:MatrixZonotopeProduct{NM}}
    return _overapproximate_lmzp(lm)
end

function _compute_inner_powers(MZ::MatrixZonotope, P::S,
                               k::Int) where {S<:Union{SparsePolynomialZonotope,
                                                       AbstractZonotope}}
    N = promote_type(eltype(MZ), eltype(P))
    S_type = Base.typename(S).wrapper

    powers = Vector{S_type{N}}(undef, k + 1)
    invfact = N(1)
    powers[1] = P

    @inbounds for i in 1:k
        invfact /= i
        term = overapproximate(invfact * MZ * powers[i], S_type)
        powers[i + 1] = remove_redundant_generators(term)
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
            term = remove_redundant_generators(overapproximate(MZ * term, S))
        end
        out_powers[i] = term
    end
    return out_powers
end

"""
    _taylor_expmap(MZ::MatrixZonotope, P::S, k::Int) where S

Compute the k-th order truncated Taylor expansion of the the exponential mao of a matrix zonotope
with a zonotopic set following Proposition 3 of [HuangLBS2025](@citet).

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
function _taylor_expmap(MZ::MatrixZonotope, P::S,
                        k::Int) where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}
    # avoid seg fault
    if isempty(generators(MZ))
        σ = linear_map(exp(center(MZ)), P)
        return remove_redundant_generators(σ)
    end

    inner = _compute_inner_powers(MZ, P, k)
    σ = reduce(exact_sum, inner)
    return remove_redundant_generators(σ)
end

"""
    _taylor_expmap(MZP::MatrixZonotopeProduct, P::S, k::Int) where S

Compute the k-th order truncated Taylor expansion of the the exponential map of the product of
matrix zonotopes with a zonotopic set following Proposition 3 of [HuangLBS2025](@citet).

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
function _taylor_expmap(MZP::MatrixZonotopeProduct, P::S,
                        k::Int) where {S<:Union{SparsePolynomialZonotope,AbstractZonotope}}
    # inner powers on the last factor
    MZP = remove_redundant_factors(MZP)
    last_factor = factors(MZP)[end]
    terms = _compute_inner_powers(last_factor, P, k)

    # propagate outwards through the preceding factors
    for MZ in reverse(view(factors(MZP), 1:(nfactors(MZP) - 1)))
        terms = _compute_outer_powers(MZ, terms, k)
    end

    σ = reduce(exact_sum, terms)
    return remove_redundant_generators(σ)
end

function load_intervalmatrices_overapproximation_matrixzonotope()
    return quote
        using .IntervalMatrices: IntervalMatrix

        function _overapproximate_emz(em::ExponentialMap{N,S,NM,MAT},
                                  k::Int=2) where {N,
                                                   S<:Union{SparsePolynomialZonotope{N},
                                                            AbstractZonotope},
                                                   NM,
                                                   MAT<:AbstractMatrixZonotope{NM}}
            T = promote_type(N, NM)
            MZP = matrix(em)
            P = set(em) #SPZ or Zonotope
            n = size(MZP, 2)

            matnorm = T(overapproximate_norm(MZP, Inf))
            ϵ = matnorm / (k + 2)
            if ϵ > 1
                @warn "k should be chosen such that ϵ<1 " ϵ
            end

            σ = _taylor_expmap(MZP, P, k)
            ε = IntervalMatrix(fill(IA.interval(T(-1), T(1)), n, n))
            ε *= matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

            Zp = overapproximate(P, Zonotope)
            rhs = overapproximate(ε * Zp, Zonotope)
            P_approx = minkowski_sum(σ, rhs)

            return remove_redundant_generators(P_approx)
        end

        """
        	function overapproximate(em::ExponentialMap{N,S,NM,MAT},
                                 T::Type{<:SparsePolynomialZonotope},
                                 k::Int=2) where {N,S<:SparsePolynomialZonotope,NM,
                                                  MAT<:AbstractMatrixZonotope}

        Overapproximate the exponential map of a sparse polynomial zonotope through a composition of matrix 
        zonotopes, following Proposition 1 of [HuangLBS2025](@citet).

        ### Input

        - `em` -- an expontial map of a sparse polynomial zonotope through a product of matrix zonotopes

        ### Output

        A sparse polynomial zonotope overapproximating the exponential map.
        """
        function overapproximate(em::ExponentialMap{N,S,NM,MAT},
                                 T::Type{<:SparsePolynomialZonotope},
                                 k::Int=2) where {N,S<:SparsePolynomialZonotope,NM,
                                                  MAT<:AbstractMatrixZonotope}
            return _overapproximate_emz(em, k)
        end

        """
        	function overapproximate(em::ExponentialMap{N,S,NM,MAT},
                                 T::Type{<:Zonotope},
                                 k::Int=2) where {N,S<:AbstractZonotope,NM,
                                                  MAT<:AbstractMatrixZonotope}

        Overapproximate the exponential map of a zonotope through a composition of matrix 
        zonotopes, following Proposition 1 of [HuangLBS2025](@citet).

        ### Input

        - `em` -- an expontial map of a zonotope through a product of matrix zonotopes

        ### Output

        A zonotope overapproximating the exponential map.
        """
        function overapproximate(em::ExponentialMap{N,S,NM,MAT}, T::Type{<:Zonotope},
                                 k::Int=2) where {N,S<:AbstractZonotope,NM,
                                                  MAT<:AbstractMatrixZonotope}
            return _overapproximate_emz(em, k)
        end
    end
end
