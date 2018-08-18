export AbstractPolytope,
       vertices_list,
       singleton_list,
       linear_map

"""
    AbstractPolytope{N<:Real} <: LazySet{N}

Abstract type for polytopic sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces,
or equivalently, sets with finitely many vertices.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

```jldoctest
julia> subtypes(AbstractPolytope)
4-element Array{Any,1}:
 AbstractPointSymmetricPolytope
 AbstractPolygon
 HPolytope
 VPolytope
```
"""
abstract type AbstractPolytope{N<:Real} <: LazySet{N} end


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
    return [Singleton(vi) for vi in vertices_list(P)]
end

"""
    linear_map(M::AbstractMatrix, P::AbstractPolytope{N}) where {N<:Real}

Concrete linear map of an abstract polytype.

### Input

- `M` -- matrix
- `P` -- abstract polytype

### Output

The polytope in V-representation obtained by applying the linear map ``M`` to
the set ``P``. If the given polytope is two-dimensional, a polygon instead
of a general polytope is returned. 
"""
function linear_map(M::AbstractMatrix, P::AbstractPolytope{N}) where {N<:Real}
    @assert dim(P) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(P))"

    if dim(P) == 2
        T = VPolygon
    else
        T = VPolytope
    end
    vlist = vertices_list(P)
    new_vlist = Vector{Vector{N}}(undef, length(vlist))
    @inbounds for (i, vi) in enumerate(vlist)
        new_vlist[i] =  M * vi
    end
    return T(new_vlist)
end
