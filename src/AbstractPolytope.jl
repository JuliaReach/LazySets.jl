import Base.isempty

export AbstractPolytope,
       vertices_list,
       singleton_list,
       isempty

"""
    AbstractPolytope{N<:Real} <: AbstractPolyhedron{N}

Abstract type for polytopic sets, i.e., bounded sets with finitely many flat
facets, or equivalently, bounded sets defined as an intersection of a finite
number of half-spaces, or equivalently, bounded sets with finitely many
vertices.

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
abstract type AbstractPolytope{N<:Real} <: AbstractPolyhedron{N} end


# =============================================
# Common AbstractPolytope functions
# =============================================

"""
    isbounded(P::AbstractPolytope)::Bool

Determine whether a polytopic set is bounded.

### Input

- `P` -- polytopic set

### Output

`true` (since a polytope must be bounded).
"""
function isbounded(::AbstractPolytope)::Bool
    return true
end

"""
    singleton_list(P::AbstractPolytope{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polytopic set as a list of singletons.

### Input

- `P` -- polytopic set

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::AbstractPolytope{N}
                       )::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in vertices_list(P)]
end

"""
    isempty(P::AbstractPolytope)::Bool

Determine whether a polytope is empty.

### Input

- `P` -- abstract polytope

### Output

`true` if the given polytope contains no vertices, and `false` otherwise.

### Algorithm

This algorithm checks whether the `vertices_list` of the given polytope is empty
or not.
"""
function isempty(P::AbstractPolytope)::Bool
    return isempty(vertices_list(P))
end

function default_polyhedra_backend(P, N)
    @assert isdefined(@__MODULE__, :Polyhedra) "this function needs the package 'Polyhedra' to be loaded"
    error("no default backend for numeric type $N")
end

# given a polytope P, apply the linear map P to each vertex of P
# it is assumed that the interface function `vertices_list(P)` is available 
@inline function _linear_map_vrep(M::AbstractMatrix{N}, P::AbstractPolytope{N}) where {N<:Real}
    vertices = broadcast(v -> M * v, vertices_list(P))
    m = size(M, 1) # output dimension
    if m == 1
        # TODO: substitute with Interval(convex_hull(vertices)...): convex hull 1D
        return Interval(minimum(vertices)[1], maximum(vertices)[1])
    elseif m == 2
        return VPolygon(vertices)
    else
        return VPolytope(vertices)
    end
end

@inline function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolytope{N},
                          use_inv::Bool) where {N<:Real}
    constraints = _linear_map_hrep_helper(M, P, use_inv)
    m = size(M, 1) # output dimension
    if m == 1
        return convert(Interval, HPolygon(constraints))
    elseif m == 2
        return HPolygon(constraints)
    else
        return HPolytope(constraints)
    end
end

# =============================================
# Functions that require Polyhedra
# =============================================

function load_polyhedra_abstractpolytope() # function to be loaded by Requires

return quote

@static if VERSION < v"0.7-"
    import Polyhedra:polyhedron,
                     HRep, VRep,
                     removehredundancy!, removevredundancy!,
                     hreps, vreps,
                     intersect,
                     convexhull,
                     hcartesianproduct, vcartesianproduct,
                     points

    import CDDLib # default backend

    function default_polyhedra_backend(P, N::Type{<:AbstractFloat})
        return CDDLib.CDDLibrary()
    end

    function default_polyhedra_backend(P, N::Type{<:Rational})
        return CDDLib.CDDLibrary(:exact)
    end

else
    import .Polyhedra:polyhedron,
                     HRep, VRep,
                     removehredundancy!, removevredundancy!,
                     hreps, vreps,
                     intersect,
                     convexhull,
                     hcartesianproduct, vcartesianproduct,
                     points,
                     default_library

    function default_polyhedra_backend(P, N::Type{<:AbstractFloat})
        return default_library(LazySets.dim(P), Float64)
    end

    function default_polyhedra_backend(P, N::Type{<:Rational})
        return default_library(LazySets.dim(P), Rational{Int})
    end
end
end # quote
end # function load_polyhedra_hpolytope()
