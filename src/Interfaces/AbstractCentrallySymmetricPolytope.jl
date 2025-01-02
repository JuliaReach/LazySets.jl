export AbstractCentrallySymmetricPolytope

"""
    AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N}

Abstract type for centrally symmetric, polytopic sets.
It combines the `AbstractCentrallySymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractCentrallySymmetricPolytope` must define the following
`AbstractCentrallySymmetric` method, in addition to the `AbstractPolytope`
methods:

- `center(::AbstractCentrallySymmetricPolytope)` -- return the center point

The subtypes of `AbstractCentrallySymmetricPolytope` (including abstract
interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetricPolytope)
2-element Vector{Any}:
 AbstractZonotope
 Ball1
```
"""
abstract type AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N} end

# common AbstractCentrallySymmetric functions
# copy-pasted because Julia does not have multiple inheritance

@inline function dim(P::AbstractCentrallySymmetricPolytope)
    return length(center(P))
end

function an_element(P::AbstractCentrallySymmetricPolytope)
    return center(P)
end

function isempty(::AbstractCentrallySymmetricPolytope)
    return false
end

function isuniversal(S::AbstractCentrallySymmetricPolytope, witness::Bool=false)
    if witness
        N = eltype(S)
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end

@inline function center(S::AbstractCentrallySymmetricPolytope, i::Int)
    return center(S)[i]
end

function extrema(S::AbstractCentrallySymmetricPolytope)
    # h = c + r
    h = high(S)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 .* center(S) .- h
    return (l, h)
end

function extrema(S::AbstractCentrallySymmetricPolytope, i::Int)
    # h = c + r
    h = high(S, i)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 * center(S, i) - h
    return (l, h)
end
