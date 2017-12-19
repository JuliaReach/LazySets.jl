export AbstractPointSymmetric,
       center

"""
    AbstractPointSymmetric{N<:Real} <: LazySet

Abstract type for point symmetric sets.

### Notes

Every concrete `AbstractPointSymmetric` must define the following functions:
- `center(::AbstractPointSymmetric{N})::Vector{N}` -- return the center point

```jldoctest
julia> subtypes(AbstractPointSymmetric)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.Ball2                   
 LazySets.Ballp                   
```
"""
abstract type AbstractPointSymmetric{N<:Real} <: LazySet end


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
