export SparsePolynomialZonotope, expmat, nparams, ngens_dep, ngens_indep,
       genmat_dep, genmat_indep, indexvector, exact_sum, ⊞, translate,
       linear_map, quadratic_map, remove_redundant_generators, reduce_order

"""
    SparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, MNI<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}, VI<:AbstractVector{<:Integer}} <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope.

A sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is represented by the set
```math
\\mathcal{PZ} = \\left\\{x \\in \\mathbb{R}^n : x = c + ∑ᵢ₌₁ʰ\\left(∏ₖ₌₁ᵖ α_k^{E_{k, i}} \\right)Gᵢ+∑ⱼ₌₁^qβⱼGIⱼ,~~ α_k ∈ [-1, 1]~~ ∀ i = 1,…,p, j=1,…,q \\right\\},
```
where ``c ∈ \\mathbb{R}^n`` is the offset vector (or center), ``G ∈ \\mathbb{R}^{n \\times h}`` is the dependent generator matrix with columns ``Gᵢ``,
``GI ∈ \\mathbb{R}^{n×q}`` is the independent generator matrix and ``E ∈ \\mathbb{N}^{p×h}_{≥0}`` is the exponent matrix with matrix elements ``E_{k, i}``.

### Fields

- `c`   -- offset vector
- `G`   -- dependent generator matrix
- `GI`  -- independent generator matrix
- `E`   -- exponent matrix
- `idx` -- identifier vector, vector of positive integers identifing the dependent parameters of `PZ`.

### Notes

Sparse polynomial zonotopes were introduced in [KA21].

- [KA21] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*. Transactions on Automatic Control, 2021.
"""
struct SparsePolynomialZonotope{N,
                                VN<:AbstractVector{N},
                                MN<:AbstractMatrix{N},
                                MNI<:AbstractMatrix{N},
                                ME<:AbstractMatrix{<:Integer},
                                VI<:AbstractVector{<:Integer}} <: AbstractPolynomialZonotope{N}
    c::VN
    G::MN
    GI::MNI
    E::ME
    idx::VI

    function SparsePolynomialZonotope(c::VN, G::MN, GI::MNI, E::ME, idx::VI) where {N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, MNI<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}, VI<:AbstractVector{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G should have the same number of rows"))
        @assert length(c) == size(GI, 1) throw(DimensionMismatch("c and GI should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should have non-negative integers"))
        @assert all(>(0), idx) throw(ArgumentError("identifiers in index vector must be positive integers"))
        return new{N, VN, MN, MNI, ME, VI}(c, G, GI, E, idx)
    end
end

function SparsePolynomialZonotope(c::AbstractVector, G::AbstractMatrix, GI::AbstractMatrix, E::AbstractMatrix{<:Integer})
    n = size(E, 1)
    return SparsePolynomialZonotope(c, G, GI, E, uniqueID(n))
end


const SPZ = SparsePolynomialZonotope

"""
    dim(P::SparsePolynomialZonotope)

Return the dimension of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The ambient dimension of `P`.
"""
dim(P::SPZ) = length(P.c)

"""
    ngens_dep(P::SparsePolynomialZonotope)

Return the number of dependent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
ngens_dep(P::SPZ) = size(P.G, 2)

"""
    ngens_indep(P::SparsePolynomialZonotope)

Return the number of independent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
ngens_indep(P::SPZ) = size(P.GI, 2)

"""
    nparams(P::SparsePolynomialZonotope)

Return the number of dependent parameters in the polynomial representation of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of dependent parameters in the polynomial representation of P.

### Notes

This corresponds to the number rows in the exponent matrix ``E``.
"""
nparams(P::SPZ) = size(P.E, 1)

"""
    order(P::SparsePolynomialZonotope)

Return the order of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The order of `P`, defined as the quotient between the number of generators and the ambient dimension.
"""
order(P::SPZ) = (ngens_dep(P) + ngens_indep(P)) // dim(P)

"""
    center(P::SparsePolynomialZonotope)

Return the center of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The center of `P`.
"""
center(P::SPZ) = P.c

"""
    genmat_dep(P::SparsePolynomialZonotope)

Return the matrix of dependent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
genmat_dep(P::SPZ) = P.G

"""
    genmat_indep(P::SparsePolynomialZonotope)

Return the matrix of independent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
genmat_indep(P::SPZ) = P.GI

"""
    expmat(P::SparsePolynomialZonotope)

Return the matrix of exponents of the sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``\alpha_k`` in the definition) and each column to a monomial.
"""
expmat(P::SPZ) = P.E

"""
    indexvector(P::SparsePolynomialZonotope)

Return the index vector of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Notes

The index vector is a vector of positive integers identifing the dependent parameters of `P`.
"""
indexvector(P::SPZ) = P.idx


"""
    uniqueID(n::Int)

Returns a collection of n unique identifiers (intergers 1, …, n).
"""
uniqueID(n::Int) = 1:n
"""
    linear_map(M::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling}, P::SparsePolynomialZonotope)

Apply the linear map `M` to the sparse polynomial zonotope `P`.

### Input

- `M` -- square matrix with `size(M) == dim(P)`
- `P` -- sparse polynomial zonotope

### Output

The set resulting from applying the linear map `M` to `P`.
"""
function linear_map(M::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling}, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * genmat_dep(P),
                                    M * genmat_indep(P),
                                    expmat(P),
                                    indexvector(P)
                                    )
end

"""
    rand(::Type{SparsePolynomialZonotope};
         [N]::Type{<:Real}=Float64, [dim]::Int=2, [nparams]::Int=2, [maxdeg]::Int=3,
         [num_dependent_generators]::Int=-1,
         [num_independent_generators]::Int=-1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random sparse polynomial zonotope.

### Input

- `Zonotope`                   -- type for dispatch
- `N`                          -- (optional, default: `Float64`) numeric type
- `dim`                        -- (optional, default: 2) dimension
- `nparams`                    -- (optional, default: 2) number of parameters
- `maxdeg`                     -- (optinal, default: 3) maximum degree for each parameter
- `rng`                        -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`                       -- (optional, default: `nothing`) seed for reseeding
- `num_dependent_generators`   -- (optional, default: `-1`) number of dependent generators of
                                  the zonotope (see comment below)
- `num_independent_generators` -- (optional, default: `-1`) number of independent generators of
                                  the zonotope (see comment below)

### Output

A random sparse polynomial zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the arguments `num_dependent_generators` and `num_dependent_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator). Note that the final number of
generators may be lower if redundant monomials are generated.
"""
function rand(::Type{SparsePolynomialZonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              nparams::Int=2,
              maxdeg::Int=3,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_dependent_generators::Int=-1,
              num_independent_generators::Int=-1)
    rng = reseed(rng, seed)

    if num_independent_generators < 0
        num_independent_generators = (dim == 1) ? 1 : rand(dim:2*dim)
    end
    GI = randn(rng, N, dim, num_independent_generators)

    SSPZ = rand(SimpleSparsePolynomialZonotope; N=N, dim=dim, nparams=nparams, maxdeg=maxdeg, rng=rng, seed=seed, num_generators=num_dependent_generators)

    return SparsePolynomialZonotope(center(SSPZ), genmat(SSPZ), GI, expmat(SSPZ))
end


"""
    exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

Compute the exact sum of sparse polyomial zonotopes ``P₁`` and ``P₂``.

### Input

- `P1` -- sparse polynomial zonotope
- `P2` -- sparse polynomial zonotope

### Output

exact sum ``P₁ ⊞ P₂``

"""
function exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

    indexvector(P1) == indexvector(P2) || throw(ArgumentError("Exact sum currently only implemented for Sparse Polynomial Zonotopes with the same index vector"))

    c = center(P1) + center(P2)
    G = hcat(genmat_dep(P1), genmat_dep(P2))
    GI = hcat(genmat_indep(P1), genmat_indep(P2))
    E = hcat(expmat(P1), expmat(P2))
    idx = indexvector(P1)

    return SparsePolynomialZonotope(c, G, GI, E, idx)
end

const ⊞ = exact_sum

"""
    translate(S::SparsePolynomialZonotope, v::AbstractVector)

Translate (i.e., shift) a sparse polynomial zonotope by a given vector.

### Input

- `S` -- sparse polynomial zonotope
- `v` -- translation vector

### Output

A translated sparse polynomial zonotope.
"""
function translate(S::SparsePolynomialZonotope, v::AbstractVector)
    c = center(S) + v
    return SparsePolynomialZonotope(c, genmat_dep(S), genmat_indep(S), expmat(S))
end

"""
    remove_redundant_generators(S::SparsePolynomialZonotope)

Remove redundant generators from `S`.

### Input

- `S` -- sparse polynomial zonotope

### Output

A new sparse polynomial zonotope such that redundant generators have been removed.

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the dependent generator matrix, `E` the exponent matrix and `GI` the independent generator matrix of `S`.
The following simplifications are performed:

- Zero columns in `G` and the corresponding columns in `E` are removed.
- Zero columns in `GI` are removed.
- For zero columns in `E`, the corresponding column in `G` is summed to the center.
- Repeated columns in `E` are grouped together by summing the corresponding columns in `G`.
"""
function remove_redundant_generators(S::SparsePolynomialZonotope)
    c, G, E = _remove_redundant_generators_polyzono(center(S), genmat_dep(S), expmat(S))
    GI = remove_zero_columns(genmat_indep(S))
    return SparsePolynomialZonotope(c, G, GI, E)
end

"""
    reduce_order(P::SparsePolynomialZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

Overapproximate the sparse polynomial zonotope `P` by one which has at most order `r`.

### Input

- `P`       -- sparse polynomial zonotope
- `r`       -- maximum order of the resulting sparse polynomial zonotope (≥ 1)
- `method`  -- (optional default [`GIR05`](@ref)) algorithm used internally for
               the order reduction of a (normal) zonotope

### Output

A sparse polynomial zonotope with order at most `r`.

### Notes

This method implements the algorithm described in Proposition 3.1.39 of [1].

[1] N. Kochdumper. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems*. 2021.
"""
function reduce_order(P::SparsePolynomialZonotope, r::Real, method::AbstractReductionMethod=GIR05())
    @assert r ≥ 1
    n = dim(P)
    h = ngens_dep(P)
    q = ngens_indep(P)

    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    idx = indexvector(P)

    a = max(0, min(h + q, ceil(Int, h + q - n * (r - 1))))
    Gbar = hcat(G, GI)
    norms = [norm(g) for g in eachcol(Gbar)]
    th = sort(norms)[a]

    # TODO: case a = 0
    K = [norms[i] ≤ th for i in 1:h] # ? Is constructing an array of booleans the most efficient way?
    Kbar = .!K

    H = [norms[h+i] ≤ th for i in 1:q]
    Hbar = .!H

    PZ = SparsePolynomialZonotope(c, G[:, K], GI[:, H], E[:, K], idx)
    Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)

    Ebar = E[:, Kbar]
    N = [!iszero(e) for e in eachrow(Ebar)]

    cz = center(Z)
    Gz = genmat(Z)
    return SparsePolynomialZonotope(cz, G[:, Kbar], hcat(GI[:, Hbar], Gz), Ebar[N, :], idx[N])
end

# function quadratic_map(Q::Vector{MT}, S::SparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
#     m = length(Q)
#     c = center(S)
#     h = ngens(S)
#     G = genmat(S)
#     E = expmat(S)

#     cnew = similar(c, m)
#     Gnew = similar(G, m, h^2 + h)
#     QiG = similar(Q)
#     @inbounds for (i, Qi) in enumerate(Q)
#         cnew[i] = dot(c, Qi, c)
#         Gnew[i, 1:h] = c' * (Qi + Qi') * G
#         QiG[i] = Qi * G
#     end

#     Enew = repeat(E, 1, h + 1)
#     @inbounds for i in 1:h
#         idxstart = h * i + 1
#         idxend = (i + 1) * h
#         Enew[:, idxstart:idxend] .+= E[:, i]
#         for j in eachindex(QiG)
#             Gnew[j, idxstart:idxend] = G[:, i]' * QiG[j]
#         end
#     end
#     return SparsePolynomialZonotope(cnew, Gnew, Enew)
# end
