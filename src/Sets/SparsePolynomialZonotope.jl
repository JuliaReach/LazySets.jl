export SparsePolynomialZonotope, expmat, nparams, ngens_dep, ngens_indep,
       genmat_dep, genmat_indep, indexvector, translate,
       linear_map, quadratic_map, remove_redundant_generators, reduce_order, ρ

"""
    SparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N},
                             MNI<:AbstractMatrix{N},
                             ME<:AbstractMatrix{<:Integer},
                             VI<:AbstractVector{<:Integer}}
        <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope.

A sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is represented by the set
```math
\\mathcal{PZ} = \\left\\{x \\in \\mathbb{R}^n : x = c + ∑ᵢ₌₁ʰ\\left(∏ₖ₌₁ᵖ α_k^{E_{k, i}} \\right)Gᵢ+∑ⱼ₌₁^qβⱼGIⱼ,~~ α_k, βⱼ ∈ [-1, 1],~~ ∀ k = 1,…,p, j=1,…,q \\right\\},
```
where ``c ∈ \\mathbb{R}^n`` is the offset vector (or center),
``Gᵢ ∈ \\mathbb{R}^{n}`` are the dependent generators,
``GIⱼ ∈ \\mathbb{R}^{n}`` are the independent generators, and
``E ∈ \\mathbb{N}^{p×h}_{≥0}`` is the exponent matrix with matrix elements ``E_{k, i}``.

In the implementation, ``Gᵢ ∈ \\mathbb{R}^n`` are arranged as columns of the dependent generator
matrix ``G ∈ \\mathbb{R}^{n \\times h}``, and similarly ``GIⱼ ∈ \\mathbb{R}^{n}`` are arranged as
columns of the independent generator matrix ``GI ∈ \\mathbb{R}^{n×q}``.

The shorthand notation ``\\mathcal{PZ} = \\langle c, G, GI, E, idx \\rangle`` is often used, where
``idx ∈ \\mathbb{N}^p`` is a list of non-repeated natural numbers
storing a unique identifier for each dependent factor ``αₖ``.

### Fields

- `c`   -- offset vector
- `G`   -- dependent generator matrix
- `GI`  -- independent generator matrix
- `E`   -- exponent matrix
- `idx` -- identifier vector of positive integers for the dependent parameters

### Notes

Sparse polynomial zonotopes were introduced in [1].

- [1] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*.
      Transactions on Automatic Control, 2021.
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

    # default constructor with dimension checks
    function SparsePolynomialZonotope(c::VN, G::MN, GI::MNI, E::ME,
                                      idx::VI=uniqueID(size(E, 1))) where {N,VN<:AbstractVector{N},
                                                                           MN<:AbstractMatrix{N},
                                                                           MNI<:AbstractMatrix{N},
                                                                           ME<:AbstractMatrix{<:Integer},
                                                                           VI<:AbstractVector{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G " *
                                                                "should have the same number of rows"))
        @assert length(c) == size(GI, 1) throw(DimensionMismatch("c and GI " *
                                                                 "should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E " *
                                                                 "should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should contain " *
                                                  "non-negative integers"))
        @assert all(>(0), idx) throw(ArgumentError("identifiers in index " *
                                                   "vector must be positive integers"))
        return new{N,VN,MN,MNI,ME,VI}(c, G, GI, E, idx)
    end
end

# short-hand
const SPZ = SparsePolynomialZonotope

function isoperationtype(P::Type{<:SparsePolynomialZonotope})
    return false
end

"""
    ngens_dep(P::SparsePolynomialZonotope)

Return the number of dependent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of dependent generators.
"""
ngens_dep(P::SPZ) = size(P.G, 2)

"""
    ngens_indep(P::SparsePolynomialZonotope)

Return the number of independent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of independent generators.
"""
ngens_indep(P::SPZ) = size(P.GI, 2)

"""
    nparams(P::SparsePolynomialZonotope)

Return the number of dependent parameters in the polynomial representation of a
sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of dependent parameters in the polynomial representation.

### Notes

This number corresponds to the number of rows in the exponent matrix ``E``.
"""
nparams(P::SPZ) = size(P.E, 1)

"""
    order(P::SparsePolynomialZonotope)

Return the order of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The order, defined as the quotient between the number of generators and the
ambient dimension, as a `Rational` number.
"""
order(P::SPZ) = (ngens_dep(P) + ngens_indep(P)) // dim(P)

"""
    center(P::SparsePolynomialZonotope)

Return the center of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The center.
"""
center(P::SPZ) = P.c

"""
    genmat_dep(P::SparsePolynomialZonotope)

Return the matrix of dependent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of dependent generators.
"""
genmat_dep(P::SPZ) = P.G

"""
    genmat_indep(P::SparsePolynomialZonotope)

Return the matrix of independent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of independent generators.
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

In the exponent matrix, each row corresponds to a parameter (``αₖ`` in the
definition) and each column to a monomial.
"""
expmat(P::SPZ) = P.E

"""
    indexvector(P::SparsePolynomialZonotope)

Return the index vector of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The index vector.

### Notes

The index vector contains positive integers for the dependent parameters.
"""
indexvector(P::SPZ) = P.idx

"""
    uniqueID(n::Int)

Return a collection of n unique identifiers (integers 1, …, n).

### Input

- `n` -- number of variables

### Output

`1:n`.
"""
uniqueID(n::Int) = 1:n

"""
    linear_map(M::AbstractMatrix, P::SparsePolynomialZonotope)

Apply a linear map to a sparse polynomial zonotope.

### Input

- `M` -- square matrix with `size(M) == dim(P)`
- `P` -- sparse polynomial zonotope

### Output

The sparse polynomial zonotope resulting from applying the linear map.
"""
function linear_map(M::AbstractMatrix, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * genmat_dep(P),
                                    M * genmat_indep(P),
                                    expmat(P),
                                    indexvector(P))
end

"""
    rand(::Type{SparsePolynomialZonotope}; [N]::Type{<:Real}=Float64,
         [dim]::Int=2, [nparams]::Int=2, [maxdeg]::Int=3,
         [num_dependent_generators]::Int=-1,
         [num_independent_generators]::Int=-1, [rng]::AbstractRNG=GLOBAL_RNG,
         [seed]::Union{Int, Nothing}=nothing)

Create a random sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope`   -- type for dispatch
- `N`                          -- (optional, default: `Float64`) numeric type
- `dim`                        -- (optional, default: 2) dimension
- `nparams`                    -- (optional, default: 2) number of parameters
- `maxdeg`                     -- (optional, default: 3) maximum degree for each
                                  parameter
- `num_dependent_generators`   -- (optional, default: `-1`) number of dependent
                                  generators (see comment below)
- `num_independent_generators` -- (optional, default: `-1`) number of
                                  independent generators (see comment below)
- `rng`                        -- (optional, default: `GLOBAL_RNG`) random
                                  number generator
- `seed`                       -- (optional, default: `nothing`) seed for
                                  reseeding

### Output

A random sparse polynomial zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the arguments
`num_dependent_generators` and `num_dependent_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator). Note that the
final number of generators may be lower if redundant monomials are generated.
"""
function rand(::Type{SparsePolynomialZonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              nparams::Int=2,
              maxdeg::Int=3,
              num_dependent_generators::Int=-1,
              num_independent_generators::Int=-1,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)

    if num_independent_generators < 0
        num_independent_generators = (dim == 1) ? 1 : rand(rng, dim:(2 * dim))
    end
    GI = randn(rng, N, dim, num_independent_generators)

    SSPZ = rand(SimpleSparsePolynomialZonotope; N=N, dim=dim, nparams=nparams,
                maxdeg=maxdeg, rng=rng, num_generators=num_dependent_generators)

    return SparsePolynomialZonotope(center(SSPZ), genmat(SSPZ), GI, expmat(SSPZ))
end

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

A new sparse polynomial zonotope where redundant generators have been removed.

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the dependent generator matrix, `E` the exponent matrix, and `GI` the
independent generator matrix of `S`. We perform the following simplifications:

- Remove zero columns in `G` and the corresponding columns in `E`.
- Remove Zero columns in `GI`.
- For zero columns in `E`, add the corresponding column in `G` to the center.
- Group repeated columns in `E` together by summing the corresponding columns in
  `G`.
"""
function remove_redundant_generators(S::SparsePolynomialZonotope)
    c, G, E = _remove_redundant_generators_polyzono(center(S), genmat_dep(S),
                                                    expmat(S))
    GI = remove_zero_columns(genmat_indep(S))
    return SparsePolynomialZonotope(c, G, GI, E)
end

"""
    reduce_order(P::SparsePolynomialZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

Overapproximate the sparse polynomial zonotope by another sparse polynomial
zonotope with order at most `r`.

### Input

- `P`       -- sparse polynomial zonotope
- `r`       -- maximum order of the resulting sparse polynomial zonotope (≥ 1)
- `method`  -- (optional default [`GIR05`](@ref)) algorithm used internally for
               the order reduction of a (normal) zonotope

### Output

A sparse polynomial zonotope with order at most `r`.

### Notes

This method implements the algorithm described in Proposition 3.1.39 of [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application to verification of cyber-physical systems.*
    PhD diss., Technische Universität München, 2022.
"""
function reduce_order(P::SparsePolynomialZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    @assert r ≥ 1 "cannot reduce below order 1 (got $r)"

    if order(P) <= r
        return P
    end

    n = dim(P)
    h = ngens_dep(P)
    q = ngens_indep(P)

    a = min(h + q, ceil(Int, h + q - n * (r - 1)))
    @assert a > 0  # holds because `r > order(P)`

    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    idx = indexvector(P)

    Gbar = hcat(G, GI)
    norms = [norm(g) for g in eachcol(Gbar)]
    th = sort(norms)[a]

    # TODO is constructing an array of booleans the most efficient way?
    K = [norms[i] ≤ th for i in 1:h]
    Kbar = .!K

    H = [norms[h + i] ≤ th for i in 1:q]
    Hbar = .!H

    PZ = SparsePolynomialZonotope(c, G[:, K], GI[:, H], E[:, K], idx)
    Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)

    Ebar = E[:, Kbar]
    N = [!iszero(e) for e in eachrow(Ebar)]

    cz = center(Z)
    Gz = genmat(Z)
    return SparsePolynomialZonotope(cz, G[:, Kbar], hcat(GI[:, Hbar], Gz),
                                    Ebar[N, :], idx[N])
end

"""
    ρ(d::AbstractVector, P::SparsePolynomialZonotope; [enclosure_method]=nothing)

Bound the support function of ``P`` in the direction ``d``.

### Input

- `d`                -- direction
- `P`                -- sparse polynomial zonotope
- `enclosure_method` -- (optional; default: `nothing`) method to use for
                        enclosure; an `AbstractEnclosureAlgorithm` from the
                        [`Rangeenclosures.jl`](https://github.com/JuliaReach/RangeEnclosures.jl)
                        package

### Output

An overapproximation of the support function in the given direction.

### Algorithm

This method implements Proposition 3.1.16 in [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application to verification of cyber-physical systems.*
    PhD diss., Technische Universität München, 2022.
"""
function ρ(d::AbstractVector, P::SparsePolynomialZonotope;
           enclosure_method=nothing)
    require(@__MODULE__, :RangeEnclosures; fun_name="ρ")
    return _ρ_range_enclosures(d, P, enclosure_method)
end

function _load_rho_range_enclosures()
    return quote
        function _ρ_range_enclosures(d::AbstractVector, P::SparsePolynomialZonotope,
                                     method::Union{RangeEnclosures.AbstractEnclosureAlgorithm,
                                                   Nothing})
            # default method: BranchAndBoundEnclosure
            isnothing(method) && (method = RangeEnclosures.BranchAndBoundEnclosure())

            c = center(P)
            G = genmat_dep(P)
            GI = genmat_indep(P)
            E = expmat(P)
            n = dim(P)

            res = d' * c + sum(abs.(d' * gi) for gi in eachcol(GI); init=zero(eltype(GI)))

            f(x) = sum(d' * gi * prod(x .^ ei) for (gi, ei) in zip(eachcol(G), eachcol(E)))

            dom = IA.IntervalBox(IA.interval(-1, 1), n)
            res += IA.sup(enclose(f, dom, method))
            return res
        end
    end
end
