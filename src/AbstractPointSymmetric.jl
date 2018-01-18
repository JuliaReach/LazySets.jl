export AbstractPointSymmetric,
       center,
       an_element

"""
    AbstractPointSymmetric{N<:Real} <: LazySet{N}

Abstract type for point symmetric sets.

### Notes

Every concrete `AbstractPointSymmetric` must define the following functions:

- `center(::AbstractPointSymmetric{N})::Vector{N}` -- return the center point

```jldoctest
julia> subtypes(AbstractPointSymmetric)
3-element Array{Union{DataType, UnionAll},1}:
 LazySets.Ball2
 LazySets.Ballp
 LazySets.Ellipsoid
```
"""
abstract type AbstractPointSymmetric{N<:Real} <: LazySet{N} end


"""
    dim(S::AbstractPointSymmetric)::Int

Return the ambient dimension of a point symmetric set.

### Input

- `S` -- set

### Output

The ambient dimension of the set.
"""
@inline function dim(S::AbstractPointSymmetric)::Int
    return length(center(S))
end

"""
    an_element(S::AbstractPointSymmetric{N})::Vector{N} where {N<:Real}

Return some element of a point symmetric set.

### Input

- `S` -- point symmetric set

### Output

The center of the point symmetric set.
"""
function an_element(S::AbstractPointSymmetric{N})::Vector{N} where {N<:Real}
    return center(S)
end
