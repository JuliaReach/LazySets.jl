"""
	overapproximate(lm::LinearMap{N, SparsePolynomialZonotope{N}, NM, MAT}) where {N, NM,
	                 MAT <: MatrixZonotope{NM}}

Overapproximate the linear map of a sparse polynomial zonotope through a matrix zonotope,
following Proposition 1 of [HuangLBS2025](@citet).

### Input

- `lm` -- a linear map of a sparse polynomial zonotope through a matrix zonotope

### Output

A sparse polynomial zonotope overapproximating the linear map.

"""
function overapproximate(lm::LinearMap{N,S,NM,MAT}) where {N,S<:SparsePolynomialZonotope{N},NM,
                                                           MAT<:MatrixZonotope{NM}}
    MZ = matrix(lm)
    P = set(lm)
    T = promote_type(N, NM)

    m, n = size(MZ)
    w = ngens(MZ)
    h = ngens_dep(P)
    q = ngens_indep(P)

    if n != dim(P)
        throw(DimensionMismatch("incompatible dimensions:" *
                                "size(MZ) = $(size(MZ)), dim(P) = $q"))
    end

    c = center(MZ) * center(P)

    # matrix of independent generators
    Gi = Matrix{T}(undef, m, q * (w + 1))
    Gi[:, 1:q] = center(MZ) * genmat_indep(P)

    # compute matrix of dependendent generators
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
    E = Matrix{eltype(Ê₁)}(undef, pₖ, h + w + h * w)
    E[:, 1:h] = Ê₁
    E[:, (h + 1):(h + w)] = Ê₂
    @inbounds for l in 1:w
        cstart = (h + w) + (l - 1) * h + 1
        cend = (h + w) + l * h
        E[:, cstart:cend] = Ê₂[:, l] * ones(1, h) .+ Ê₁
    end

    return SparsePolynomialZonotope(c, G, Gi, E, idx)
end

function overapproximate(lm::LinearMap{N,S,NM,MAT}) where {N,S<:AbstractZonotope{N},NM,
                                                           MAT<:MatrixZonotope{NM}}
    MZ = matrix(lm)
    Z = set(lm)
    T = promote_type(N, NM)

    m, n = size(MZ)
    w = ngens(MZ)
    h = ngens(Z)

    if n != dim(Z)
        throw(DimensionMismatch("incompatible dimensions:" *
                                "size(MZ) = $(size(MZ)), dim(Z) = $q"))
    end

    c = center(MZ) * center(Z)

    # compute matrix of dependendent generators
    G = Matrix{T}(undef, m, h + w + h * w)
    G[:, 1:h] = center(MZ) * genmat(Z)

    # loop to populate G 
    @inbounds for (i, A) in enumerate(generators(MZ))
        G[:, h + i] = A * center(P)
        G[:, (h + w + (i - 1) * h + 1):(h + w + i * h)] = A * genmat(Z)
    end

    return Zonotope(c, G)
end

"""
    _taylor_expmap(A::T, B::MatrixZonotope, P::SparsePolynomialZonotope, k::Int) where {T<:Union{MatrixZonotope, Nothing}}

Compute the k-th order truncated Taylor expansion of the composition of matrix zonotopes with a sparse polynomial zonotope,
following Proposition 3 of [HuangLBS2025](@citet).

### Input

- `A` -- a matrix zonotope, or `nothing` if the expansion is for a single matrix zonotope `B`
- `B` -- a matrix zonotope
- `P` -- a sparse polynomial zonotope
- `k` -- the order of the Taylor expansion

### Output

A sparse polynomial zonotope representing the k-th order truncated Taylor expansion.

### Algorithm

This function computes the approximation:

```math
\\displaystyle\\boxplus_{i=0}^k \\frac{\\mathcal{A}^i \\mathcal{B}^i}{i!} \\mathcal{PZ}
```
"""
function _taylor_expmap(A::T, B::MatrixZonotope, P::SparsePolynomialZonotope,
                        k::Int) where {T<:Union{MatrixZonotope,Nothing}}
    N = isnothing(A) ? promote_type(eltype(B), eltype(P)) :
        promote_type(eltype(A), eltype(B), eltype(P))

    #Bpowers = [P, BP, BBP/2, ..., B^k P /k!] 
    Bpowers = Vector{SparsePolynomialZonotope{N}}(undef, k + 1)
    invfact = N(1)
    Bpowers[1] = P

    @inbounds for i in 1:k
        invfact /= i
        term = overapproximate(invfact * B * Bpowers[i])
        Bpowers[i + 1] = remove_redundant_generators(term)
    end

    if isnothing(A)
        result = reduce(exact_sum, Bpowers)
        return result
    end

    # A^i * (B^i * X)/ i!
    result = Bpowers[1]
    @inbounds for i in 1:k
        term = Bpowers[i + 1]
        # apply A^i
        for _ in 1:i
            term = remove_redundant_generators(overapproximate(A * term))
        end
        result = remove_redundant_generators(exact_sum(result, term))
    end
    return result
end

# helper to pull out A,B
get_factors(MZP::MatrixZonotope) = (nothing, MZP)
get_factors(MZP::MatrixZonotopeProduct) = (MZP.A, MZP.B)

function load_intervalmatrices_overapproximation_matrixzonotope_s()
    return quote
        using .IntervalMatrices: IntervalMatrix

        """
        	overapproximate(em::ExponentialMap{N,S,NM,MAT},
        					k::Int=2) where {N,S<:SparsePolynomialZonotope{N},NM,
        									MAT<:AbstractMatrixZonotope{NM}}

        Overapproximate the exponential map of a sparse polynomial zonotope through a composition of matrix 
        zonotopes, following Proposition 1 of [HuangLBS2025](@citet).

        ### Input

        - `em` -- an expontial map of a sparse polynomial zonotope through a product of matrix zonotopes

        ### Output

        A sparse polynomial zonotope overapproximating the exponential map.
        """
        function overapproximate(em::ExponentialMap{N,S,NM,MAT},
                                 k::Int=2) where {N,S<:SparsePolynomialZonotope{N},NM,
                                                  MAT<:AbstractMatrixZonotope{NM}}
            T    = promote_type(N, NM)
            MZP  = matrix(em)
            A, B = get_factors(MZP)
            P    = set(em)
            n    = size(B, 2)

            ABn = T(A === nothing ? overapproximate_norm(B, Inf) :
                    overapproximate_norm(A, Inf) * overapproximate_norm(B, Inf)) # norm overapproximation
            ϵ = ABn / (k + 2)
            if ϵ > 1
                @warn "k should be chosen such that ϵ<1 " ϵ
            end

            σ = _taylor_expmap(A, B, P, k)
            ε = IntervalMatrix(fill(IA.interval(T(-1), T(1)), n, n))
            ε *= ABn^(k + 1) / (factorial(k + 1) * (1 - ϵ))

            Zp = overapproximate(P, Zonotope)
            rhs = overapproximate(ε * Zp, Zonotope)

            return minkowski_sum(σ, rhs)
        end
    end
end
