export AbstractPolytope,
       vertices_list,
       singleton_list

"""
    AbstractPolytope{N<:Real} <: AbstractConvexSet{N}

Abstract type for polytopic sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces,
or equivalently, sets with finitely many vertices.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

```jldoctest
julia> subtypes(AbstractPolytope)
4-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetricPolytope
 LazySets.AbstractPolygon
 LazySets.HPolytope
 LazySets.VPolytope
```
"""
abstract type AbstractPolytope{N<:Real} <: AbstractConvexSet{N} end


# --- common AbstractPolytope functions ---


"""
    singleton_list(P::AbstractPolytope{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polytopic as a list of singletons.

### Input

- `P` -- a polytopic set

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::AbstractPolytope{N}
                       )::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in P.vertices_list]
end
