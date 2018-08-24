export AbstractCentrallySymmetric,
       center,
       an_element

"""
    AbstractCentrallySymmetric{N<:Real} <: LazySet{N}

Abstract type for centrally symmetric sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following functions:

- `center(::AbstractCentrallySymmetric{N})::Vector{N}` -- return the center
    point

```jldoctest
julia> subtypes(AbstractCentrallySymmetric)
3-element Array{Any,1}:
 Ball2
 Ballp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N<:Real} <: LazySet{N} end


"""
    dim(S::AbstractCentrallySymmetric)::Int

Return the ambient dimension of a centrally symmetric set.

### Input

- `S` -- set

### Output

The ambient dimension of the set.
"""
@inline function dim(S::AbstractCentrallySymmetric)::Int
    return length(center(S))
end

"""
    an_element(S::AbstractCentrallySymmetric{N})::Vector{N} where {N<:Real}

Return some element of a centrally symmetric set.

### Input

- `S` -- centrally symmetric set

### Output

The center of the centrally symmetric set.
"""
function an_element(S::AbstractCentrallySymmetric{N})::Vector{N} where {N<:Real}
    return center(S)
end
