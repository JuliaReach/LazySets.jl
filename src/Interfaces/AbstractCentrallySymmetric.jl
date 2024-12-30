export AbstractCentrallySymmetric

"""
    AbstractCentrallySymmetric{N} <: ConvexSet{N}

Abstract type for centrally symmetric compact convex sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following function:

- `center(::AbstractCentrallySymmetric)` -- return the center point

The subtypes of `AbstractCentrallySymmetric`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetric)
2-element Vector{Any}:
 AbstractBallp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N} <: ConvexSet{N} end

@inline function dim(S::AbstractCentrallySymmetric)
    return length(center(S))
end

# a set with a unique center must be bounded
function isboundedtype(::Type{<:AbstractCentrallySymmetric})
    return true
end

# a set with a unique center must be bounded
function isbounded(::AbstractCentrallySymmetric)
    return true
end

"""
# Extended help

    an_element(S::AbstractCentrallySymmetric)

### Output

The center of the centrally symmetric set.
"""
function an_element(S::AbstractCentrallySymmetric)
    return center(S)
end

function isempty(::AbstractCentrallySymmetric)
    return false
end

"""
# Extended help

    isuniversal(S::AbstractCentrallySymmetric, [witness]::Bool=false)

### Algorithm

A witness is obtained by computing the support vector in direction
`d = [1, 0, …, 0]` and adding `d` on top.
"""
function isuniversal(S::AbstractCentrallySymmetric, witness::Bool=false)
    if witness
        N = eltype(S)
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end

@inline function center(S::AbstractCentrallySymmetric, i::Int)
    return center(S)[i]
end

"""
# Extended help

    extrema(S::AbstractCentrallySymmetric)

### Notes

The result is equivalent to `(low(S), high(S))`.

### Algorithm

We compute `high(S)` and then compute the lowest coordinates with the help of
`center(S)` (which is assumed to be cheaper to obtain).
"""
function extrema(S::AbstractCentrallySymmetric)
    # h = c + r
    h = high(S)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 .* center(S) .- h
    return (l, h)
end

"""
# Extended help

    extrema(S::AbstractCentrallySymmetric, i::Int)

### Notes

The result is equivalent to `(low(S, i), high(S, i))`.

### Algorithm

We compute `high(S, i)` and then compute the lowest coordinates with the help of
`center(S, i)` (which is assumed to be cheaper to obtain).
"""
function extrema(S::AbstractCentrallySymmetric, i::Int)
    # h = c + r
    h = high(S, i)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 * center(S, i) - h
    return (l, h)
end
