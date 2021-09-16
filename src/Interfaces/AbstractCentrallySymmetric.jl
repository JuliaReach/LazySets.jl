import Base.isempty

export AbstractCentrallySymmetric,
       center,
       an_element

"""
    AbstractCentrallySymmetric{N} <: LazySet{N}

Abstract type for centrally symmetric sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following functions:

- `center(::AbstractCentrallySymmetric)` -- return the center point
- `center(::AbstractCentrallySymmetric, i::Int)` -- return the center point at index `i`

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetric)
3-element Vector{Any}:
 Ball2
 Ballp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N} <: LazySet{N} end

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

function isboundedtype(::Type{<:AbstractCentrallySymmetric})
    return true
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
               ) where {N}

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
                    ) where {N}
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

Return the center along a given dimension of a centrally symmetric set.

### Input

- `S` -- centrally symmetric set
- `i` -- dimension of interest

### Output

The center along a given dimension of the centrally symmetric set.
"""
@inline function center(S::AbstractCentrallySymmetric, i::Int)
    return center(S)[i]
end
