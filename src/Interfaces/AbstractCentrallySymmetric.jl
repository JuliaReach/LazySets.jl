export AbstractCentrallySymmetric

"""
    AbstractCentrallySymmetric{N} <: ConvexSet{N}

Abstract type for centrally symmetric compact convex sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following functions:

- `center(::AbstractCentrallySymmetric)` -- return the center point
- `center(::AbstractCentrallySymmetric, i::Int)` -- return the center point at
                                                    index `i`

The subtypes of `AbstractCentrallySymmetric`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetric)
2-element Vector{Any}:
 AbstractBallp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N} <: ConvexSet{N} end

isconvextype(::Type{<:AbstractCentrallySymmetric}) = true

"""
    dim(S::AbstractCentrallySymmetric)

Return the ambient dimension of a centrally symmetric set.

### Input

- `S` -- centrally symmetric set

### Output

The ambient dimension of the set.
"""
@inline function dim(S::AbstractCentrallySymmetric)
    return length(center(S))
end

function isboundedtype(::Type{<:AbstractCentrallySymmetric})
    return true
end

"""
    isbounded(S::AbstractCentrallySymmetric)

Check whether a centrally symmetric set is bounded.

### Input

- `S` -- centrally symmetric set

### Output

`true` (since a set with a unique center must be bounded).
"""
function isbounded(::AbstractCentrallySymmetric)
    return true
end

"""
    an_element(S::AbstractCentrallySymmetric)

Return some element of a centrally symmetric set.

### Input

- `S` -- centrally symmetric set

### Output

The center of the centrally symmetric set.
"""
function an_element(S::AbstractCentrallySymmetric)
    return center(S)
end

"""
    isempty(S::AbstractCentrallySymmetric)

Check whether a centrally symmetric set is empty.

### Input

- `S` -- centrally symmetric set

### Output

`false`.
"""
function isempty(::AbstractCentrallySymmetric)
    return false
end

"""
    isuniversal(S::AbstractCentrallySymmetric{N},
                [witness]::Bool=false) where {N}

Check whether a centrally symmetric set is universal.

### Input

- `S`       -- centrally symmetric set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``

### Algorithm

Centrally symmetric sets are bounded.
A witness is obtained by computing the support vector in direction
`d = [1, 0, …, 0]` and adding `d` on top.
"""
function isuniversal(S::AbstractCentrallySymmetric{N},
                     witness::Bool=false) where {N}
    if witness
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end

"""
    center(H::AbstractCentrallySymmetric, i::Int)

Return the center of a centrally symmetric set along a given dimension.

### Input

- `S` -- centrally symmetric set
- `i` -- dimension of interest

### Output

The center along the given dimension.
"""
@inline function center(S::AbstractCentrallySymmetric, i::Int)
    return center(S)[i]
end

"""
    extrema(S::AbstractCentrallySymmetric)

Return two vectors with the lowest and highest coordinate of a centrally
symmetric set.

### Input

- `S` -- centrally symmetric set

### Output

Two vectors with the lowest and highest coordinates of `S`.

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
    extrema(S::AbstractCentrallySymmetric, i::Int)

Return the lower and higher coordinate of a centrally symmetric set in a given
dimension.

### Input

- `S` -- centrally symmetric set
- `i` -- dimension of interest

### Output

The lower and higher coordinate of the centrally symmetric set in the given
dimension.

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
