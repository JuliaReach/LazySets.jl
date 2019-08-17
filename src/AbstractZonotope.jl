export AbstractZonotope,
       genmat,
       generators,
       ngens

"""
    AbstractZonotope{N<:Real} <: AbstractCentrallySymmetricPolytope{N}

Abstract type for zonotopic sets.

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its center and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of generators.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``\\mathbb{R}^n`` by an affine transformation.

See [`Zonotope`](@ref) for a standard implementation of this interface.

Every concrete `AbstractZonotope` must define the following functions:
- `genmat(::AbstractZonotope{N})::AbstractMatrix{N}` -- return the generator
    matrix
- `generators(::AbstractZonotope{N})` -- return an iterator over the generators

Since the functions `genmat` and `generators` can be defined in terms of each
other, it is sufficient to only genuinely implement one of them and let the
implementation of the other function call the fallback implementation
`genmat_fallback` resp. `generators_fallback`.

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractZonotope)
3-element Array{Any,1}:
 AbstractHyperrectangle
 LineSegment
 Zonotope
```
"""
abstract type AbstractZonotope{N<:Real} <: AbstractCentrallySymmetricPolytope{N}
end


# --- fallback implementations of AbstractZonotope functions ---


"""
    genmat_fallback(Z::AbstractZonotope{N}) where {N<:Real}

Fallback definition of `genmat` for zonotopic sets.

### Input

- `Z` -- zonotopic set

### Output

A matrix where each column represents one generator of `Z`.
"""
function genmat_fallback(Z::AbstractZonotope{N}) where {N<:Real}
    gens = generators(Z)
    if isempty(gens)
        return Matrix{N}(undef, dim(Z), 0)
    end
    return hcat(gens...)
end

# iterator that wraps the generator matrix
struct FallbackGeneratorIterator{M<:AbstractMatrix}
    G::M
    n_plus_one::Int

    FallbackGeneratorIterator(G::M) where {M<:AbstractMatrix} =
        new{M}(G, size(G, 2) + 1)
end

Base.length(it::FallbackGeneratorIterator) = it.n_plus_one - 1

Base.eltype(::Type{<:FallbackGeneratorIterator{<:AbstractMatrix{N}}}) where {N} =
    AbstractVector{N}

Base.eltype(::Type{<:FallbackGeneratorIterator{<:Matrix{N}}}) where {N} =
    Vector{N}

Base.eltype(::Type{<:FallbackGeneratorIterator{<:SparseMatrixCSC{N}}}) where {N} =
    SparseVector{N}

function Base.iterate(it::FallbackGeneratorIterator, state::Int=1)
    if state == it.n_plus_one
        return nothing
    end
    g = it.G[:, state]
    state += 1
    return (g, state)
end

"""
    generators_fallback(Z::AbstractZonotope{N}) where {N<:Real}

Fallback definition of `generators` for zonotopic sets.

### Input

- `Z` -- zonotopic set

### Output

An iterator over the generators of `Z`.
"""
function generators_fallback(Z::AbstractZonotope{N}) where {N<:Real}
    return FallbackGeneratorIterator(genmat(Z))
end


# --- common AbstractZonotope functions ---


"""
    ngens(Z::AbstractZonotope)::Int

Return the number of generators of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

An integer representing the number of generators.
"""
function ngens(Z::AbstractZonotope)::Int
    return length(generators(Z))
end

"""
    minkowski_sum(Z1::AbstractZonotope{N}, Z2::AbstractZonotope{N})
        where {N<:Real}

Concrete Minkowski sum of a pair of zonotopic sets.

### Input

- `Z1` -- zonotopic set
- `Z2` -- zonotopic set

### Output

A `Zonotope` corresponding to the concrete Minkowski sum of `Z1` and `Z2`.

### Algorithm

The resulting zonotope is obtained by summing up the centers and concatenating
the generators of `Z1` and `Z2`.
"""
function minkowski_sum(Z1::AbstractZonotope{N},
                       Z2::AbstractZonotope{N}) where {N<:Real}
    return Zonotope(center(Z1) + center(Z2), [genmat(Z1) genmat(Z2)])
end

"""
    linear_map(M::AbstractMatrix{N}, Z::AbstractZonotope{N}) where {N<:Real}

Concrete linear map of a zonotopic set.

### Input

- `M` -- matrix
- `Z` -- zonotopic set

### Output

A zonotope corresponding to the concrete linear map `M` applied to `Z`.
"""
function linear_map(M::AbstractMatrix{N}, Z::AbstractZonotope{N}
                   ) where {N<:Real}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"
    return Zonotope(M * center(Z), M * genmat(Z))
end
