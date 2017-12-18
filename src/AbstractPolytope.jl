import Base.⊆

export AbstractPolytope,
       vertices_list,
       singleton_list

"""
    AbstractPolytope{N<:Real} <: LazySet

Abstract type for polytopic sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces,
or equivalently, sets with finitely many vertices.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

```jldoctest
julia> subtypes(AbstractPolytope)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetricPolytope
 LazySets.AbstractPolygon
```
"""
abstract type AbstractPolytope{N<:Real} <: LazySet end


# --- LazySet interface functions ---


"""
    ⊆(P::AbstractPolytope, S::LazySet)::Bool

Check whether a polytope is contained in a convex set.

### Input

- `P` -- polytope (containee?)
- `S` -- convex set (container?)

### Output

`true` iff ``P ⊆ S``.

### Algorithm

Since ``S`` is convex, ``P ⊆ S`` iff ``v_i ∈ S`` for all vertices ``v_i`` of
``P``.
"""
function ⊆(P::AbstractPolytope, S::LazySet)::Bool
    @assert dim(P) == dim(S)

    for v in vertices_list(P)
        if !∈(v, S)
            return false
        end
    end
    return true
end


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
