import Base.isempty

export AbstractCentrallySymmetric,
       center,
       an_element

"""
    AbstractCentrallySymmetric{N<:Real} <: LazySet{N}

Abstract type for centrally symmetric sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following functions:

- `center(::AbstractCentrallySymmetric{N})` -- return the center
    point
- `center(::AbstractCentrallySymmetric{N}, i::Int)` -- return the center point at index `i`

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetric)
3-element Array{Any,1}:
 Ball2
 Ballp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N<:Real} <: LazySet{N} end

isconvextype(::Type{<:AbstractCentrallySymmetric}) = false

"""
    dim(S::AbstractCentrallySymmetric)

Return the ambient dimension of a centrally symmetric set.

### Input

- `S` -- set

### Output

The ambient dimension of the set.
"""
@inline function dim(S::AbstractCentrallySymmetric)
    return length(center(S))
end

"""
    isbounded(S::AbstractCentrallySymmetric)

Determine whether a centrally symmetric set is bounded.

### Input

- `S` -- centrally symmetric set

### Output

`true` (since a set with a unique center must be bounded).
"""
function isbounded(::AbstractCentrallySymmetric)
    return true
end

"""
    an_element(S::AbstractCentrallySymmetric{N}) where {N<:Real}

Return some element of a centrally symmetric set.

### Input

- `S` -- centrally symmetric set

### Output

The center of the centrally symmetric set.
"""
function an_element(S::AbstractCentrallySymmetric{N}) where {N<:Real}
    return center(S)
end

"""
    isempty(S::AbstractCentrallySymmetric)

Return if a centrally symmetric set is empty or not.

### Input

- `S` -- centrally symmetric set

### Output

`false`.
"""
function isempty(::AbstractCentrallySymmetric)
    return false
end

"""
    isuniversal(S::AbstractCentrallySymmetric{N}, [witness]::Bool=false
               ) where {N<:Real}

Check whether a centrally symmetric set is universal.

### Input

- `S`       -- centrally symmetric set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``

### Algorithm

A witness is obtained by computing the support vector in direction
`d = [1, 0, …, 0]` and adding `d` on top.
"""
function isuniversal(S::AbstractCentrallySymmetric{N}, witness::Bool=false
                    ) where {N<:Real}
    if witness
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end
