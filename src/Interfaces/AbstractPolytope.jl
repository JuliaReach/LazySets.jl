import Base.isempty

export AbstractPolytope,
       vertices_list,
       isempty,
       minkowski_sum

"""
    AbstractPolytope{N<:Real} <: AbstractPolyhedron{N}

Abstract type for compact convex polytopic sets.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})` -- return a list of all vertices

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolytope)
4-element Array{Any,1}:
 AbstractCentrallySymmetricPolytope
 AbstractPolygon
 HPolytope
 VPolytope
```

A polytope is a bounded polyhedron (see [`AbstractPolyhedron`](@ref)).
Polytopes are compact convex sets with either of the following equivalent
properties:
1. They are the intersection of a finite number of closed half-spaces.
2. They are the convex hull of finitely many vertices.
"""
abstract type AbstractPolytope{N<:Real} <: AbstractPolyhedron{N} end

isconvextype(::Type{<:AbstractPolytope}) = true

# =============================================
# Common AbstractPolytope functions
# =============================================

"""
    isbounded(P::AbstractPolytope)

Determine whether a polytopic set is bounded.

### Input

- `P` -- polytopic set

### Output

`true` (since a polytope must be bounded).
"""
function isbounded(::AbstractPolytope)
    return true
end

"""
    isempty(P::AbstractPolytope)

Determine whether a polytope is empty.

### Input

- `P` -- abstract polytope

### Output

`true` if the given polytope contains no vertices, and `false` otherwise.

### Algorithm

This algorithm checks whether the `vertices_list` of the given polytope is empty
or not.
"""
function isempty(P::AbstractPolytope)
    return isempty(vertices_list(P))
end

"""
    isuniversal(P::AbstractPolytope{N}, [witness]::Bool=false
               ) where {N<:Real}

Check whether a polyhedron is universal.

### Input

- `P`       -- polyhedron
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``

### Algorithm

A witness is produced using `isuniversal(H)` where `H` is the first linear
constraint of `P`.
"""
function isuniversal(P::AbstractPolytope{N}, witness::Bool=false
                    ) where {N<:Real}
    if witness
        constraints = constraints_list(P)
        if isempty(constraints)
            return (true, N[])  # special case for polytopes without constraints
        end
        return isuniversal(constraints[1], true)
    else
        return false
    end
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

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::AbstractPolytope{N},
                                 algo::AbstractLinearMapAlgorithm) where {N<:Real}
    constraints = _linear_map_hrep(M, P, algo)
    m = size(M, 1) # output dimension
    if m == 1
        # TODO: create interval directly ?
        return convert(Interval, HPolytope(constraints))
    elseif m == 2
        return HPolygon(constraints)
    else
        return HPolytope(constraints)
    end
end

# the "backend" argument is ignored, used for dispatch
function _vertices_list(P::AbstractPolytope, backend)
    return vertices_list(P)
end
