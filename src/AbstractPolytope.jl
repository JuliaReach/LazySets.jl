import Base.isempty

export AbstractPolytope,
       vertices_list,
       singleton_list,
       isempty,
       minkowski_sum

"""
    AbstractPolytope{N<:Real} <: AbstractPolyhedron{N}

Abstract type for compact convex polytopic sets.

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

A polytope is a bounded polyhedron (see [`AbstractPolyhedron`](@ref)).
Polytopes are compact convex sets with either of the following equivalent
properties:
1. They are the intersection of a finite number of closed half-spaces.
2. They are the convex hull of finitely many vertices.
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
    require(:Polyhedra; fun_name="default_polyhedra_backend")
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

"""
    minkowski_sum(P1::AbstractPolytope{N}, P2::AbstractPolytope{N};
                  [apply_convex_hull]=true,
                  [backend]=default_polyhedra_backend(P1, N),
                  [solver]=default_lp_solver(N)) where {N<:Real}

Compute the Minkowski sum between two polytopes using their vertex representation.

### Input

- `P1`                -- polytope
- `P2`                -- another polytope
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         pairwise sums using a convex hull algorithm 
- `backend`           -- (optional, default: `default_polyhedra_backend(P1, N)`)
                         the backend for polyhedral computations used to
                         post-process with a convex hull
- `solver`            -- (optional, default: `default_lp_solver(N)`) the linear programming
                         solver used in the backend

### Output

A new polytope in vertex representation whose vertices are the convex hull of
the sum of all possible sums of vertices of `P1` and `P2`.
"""
function minkowski_sum(P1::AbstractPolytope{N}, P2::AbstractPolytope{N};
                       apply_convex_hull::Bool=true,
                       backend=nothing,
                       solver=nothing) where {N<:Real}

    @assert dim(P1) == dim(P2) "cannot compute the Minkowski sum between a polyotope " *
        "of dimension $(dim(P1)) and a polytope of dimension $((dim(P2)))"

    vlist1, vlist2 = vertices_list(P1), vertices_list(P2)
    n, m = length(vlist1), length(vlist2)
    Vout = Vector{Vector{N}}()
    sizehint!(Vout, n + m)
    for vi in vlist1
        for vj in vlist2
            push!(Vout, vi + vj)
        end
    end
    if apply_convex_hull
        if backend == nothing
            backend = default_polyhedra_backend(P1, N)
            solver = default_lp_solver(N)
        end
        convex_hull!(Vout, backend=backend, solver=solver)
    end
    return VPolytope(Vout)
end

# =============================================
# Functions that require Polyhedra
# =============================================

function load_polyhedra_abstractpolytope() # function to be loaded by Requires
return quote
import .Polyhedra: polyhedron
export polyhedron
using .Polyhedra: HRep, VRep,
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

import JuMP, GLPK

function default_lp_solver(N::Type{<:AbstractFloat})
    return JuMP.with_optimizer(GLPK.Optimizer)
end

end # quote
end # function load_polyhedra_hpolytope()
