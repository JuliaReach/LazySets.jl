export AbstractPolygon,
       tohrep,
       tovrep

"""
    AbstractPolygon{N<:Real} <: AbstractPolytope{N}

Abstract type for polygons (i.e., 2D polytopes).

### Notes

Every concrete `AbstractPolygon` must define the following functions:
- `tovrep(::AbstractPolygon{N})::VPolygon{N}`         -- transform into
    V-representation
- `tohrep(::AbstractPolygon{N})::AbstractHPolygon{N}` -- transform into
    H-representation

```jldoctest
julia> subtypes(AbstractPolygon)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractHPolygon
 LazySets.VPolygon
```
"""
abstract type AbstractPolygon{N<:Real} <: AbstractPolytope{N} end


# --- LazySet interface functions ---


"""
    dim(P::AbstractPolygon)::Int

Return the ambient dimension of a polygon.

### Input

- `P` -- polygon

### Output

The ambient dimension of the polygon, which is 2.
"""
@inline function dim(P::AbstractPolygon)::Int
    return 2
end
