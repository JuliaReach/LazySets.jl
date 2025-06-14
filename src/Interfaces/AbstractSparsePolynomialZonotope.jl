export AbstractSparsePolynomialZonotope,
       expmat,
       genmat_dep, genmat_indep,
       nparams

"""
    AbstractSparsePolynomialZonotope{N} <: AbstractPolynomialZonotope{N}

Abstract type for sparse polynomial zonotope sets.

### Notes

See [`SparsePolynomialZonotope`](@ref) for a standard implementation of this
interface.

Every concrete `AbstractSparsePolynomialZonotope` must define the following functions:

- `expmat(::AbstractSparsePolynomialZonotope)` -- return the exponent matrix (sparse PZ only)
- `genmat_dep(::AbstractSparsePolynomialZonotope)` -- return the matrix of dependent generators
- `genmat_indep(::AbstractSparsePolynomialZonotope)` -- return the matrix of independent generators

The subtypes of `AbstractSparsePolynomialZonotope` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractSparsePolynomialZonotope)
2-element Vector{Any}:
 SimpleSparsePolynomialZonotope
 SparsePolynomialZonotope
```
"""
abstract type AbstractSparsePolynomialZonotope{N} <: AbstractPolynomialZonotope{N} end

"""
    expmat(P::AbstractSparsePolynomialZonotope)

Return the matrix of exponents of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``αₖ`` in the
definition) and each column corresponds to a monomial.
"""
function expmat(::AbstractSparsePolynomialZonotope) end

"""
    genmat_dep(P::AbstractSparsePolynomialZonotope)

Return the matrix of dependent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of dependent generators.
"""
function genmat_dep(::AbstractSparsePolynomialZonotope) end

"""
    genmat_indep(P::AbstractSparsePolynomialZonotope)

Return the matrix of independent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of independent generators.
"""
function genmat_indep(::AbstractSparsePolynomialZonotope) end

function ngens_dep(P::AbstractSparsePolynomialZonotope)
    return size(genmat_dep(P), 2)
end

function ngens_indep(P::AbstractSparsePolynomialZonotope)
    return size(genmat_indep(P), 2)
end

"""
    nparams(P::AbstractSparsePolynomialZonotope)

Return the number of dependent parameters in the polynomial representation of a
sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of dependent parameters in the polynomial representation.

### Notes

This number corresponds to the number of rows in the exponent matrix.
"""
function nparams(P::AbstractSparsePolynomialZonotope)
    return size(expmat(P), 1)
end

function _remove_redundant_generators_polyzono(c, G, E)
    Gnew = Matrix{eltype(G)}(undef, size(G, 1), 0)
    Enew = Matrix{eltype(E)}(undef, size(E, 1), 0)
    cnew = copy(c)

    visited_exps = Dict{Vector{Int},Int}()
    @inbounds for (gi, ei) in zip(eachcol(G), eachcol(E))
        all(isapproxzero, gi) && continue
        if iszero(ei)
            cnew += gi
        elseif haskey(visited_exps, ei) # repeated exponent
            idx = visited_exps[ei]
            Gnew[:, idx] += gi
        else
            Gnew = hcat(Gnew, gi)
            Enew = hcat(Enew, ei)
            visited_exps[ei] = size(Enew, 2)
        end
    end

    return cnew, Gnew, Enew
end

function extrema(P::AbstractSparsePolynomialZonotope; algorithm::String="zonotope")
    if algorithm == "zonotope"
        return _extrema_polynomial_zonotope(P)
    elseif algorithm == "lowhigh"
        return _extrema_lowhigh(P)
    else
        throw(ArgumentError("unknown algorithm \"$algorithm\""))
    end
end

# See [KochdumperSAB23; Proposition 1](@citet).
function _extrema_polynomial_zonotope(P::AbstractSparsePolynomialZonotope{N}) where {N}
    G = genmat_dep(P)
    E = expmat(P)
    n = dim(P)
    g₁ = zeros(N, n)
    g₂ = zeros(N, n)
    g₃ = zeros(N, n)
    for (j, Gj) in enumerate(eachcol(G))
        if all(iseven.(@view E[:, j]))
            g₁ .+= Gj
            g₂ .+= abs.(Gj)
        else
            g₃ .+= abs.(Gj)
        end
    end
    g₄ = zeros(N, n)
    for GIj in eachcol(genmat_indep(P))
        g₄ .+= abs.(GIj)
    end
    c = center(P) + N(1 / 2) .* g₁
    r = N(1 / 2) .* g₂ + g₃ + g₄
    return (c .- r, c .+ r)
end

function linear_map(M::AbstractMatrix, P::AbstractSparsePolynomialZonotope)
    return SparsePolynomialZonotope(M * center(P),
                                    M * genmat_dep(P),
                                    M * genmat_indep(P),
                                    expmat(P))
end
