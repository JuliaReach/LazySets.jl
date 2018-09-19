import Base.isempty

export AbstractPolytope,
       vertices_list,
       singleton_list,
       linear_map,
       isempty

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
 AbstractCentrallySymmetricPolytope
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
function linear_map(M::AbstractMatrix, P::AbstractPolytope{N})::VPolytope{N} where {N<:Real}
    @assert dim(P) == size(M, 2)
    return broadcast(v -> M * v, vertices_list(P)) |> VPolytope{N}
end

"""
    isempty(P::AbstractPolytope{N})::Bool where {N<:Real}

Determine whether a polytope is empty.

### Input

- `P` -- abstract polytope

### Output

`true` if the given polytope contains no vertices, and `false` otherwise.

### Algorithm

This algorithm checks whether the `vertices_list` of the given polytope is empty
or not.
"""
function isempty(P::AbstractPolytope{N})::Bool where {N<:Real}
    return isempty(vertices_list(P))
end

function load_polyhedra_abstractpolytope() # function to be loaded by Requires

return quote

using CDDLib # default backend
import Polyhedra:polyhedron, SimpleHRepresentation, SimpleHRepresentation,
                 HRep, VRep,
                 removehredundancy!, removevredundancy!,
                 hreps, vreps,
                 intersect,
                 convexhull,
                 hcartesianproduct,
                 points

function default_polyhedra_backend(N::Type{<:AbstractFloat})
    return CDDLib.CDDLibrary()
end

function default_polyhedra_backend(N::Type{<:Rational})
    return CDDLib.CDDLibrary(:exact)
end

function default_polyhedra_backend(N)
    error("no default backend for numeric type $N")
end

end # quote
end # function load_polyhedra_hpolytope()
