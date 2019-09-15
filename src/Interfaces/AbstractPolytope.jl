import Base.isempty

export AbstractPolytope,
       vertices_list,
       singleton_list,
       isempty,
       minkowski_sum,
       chebyshev_center

"""
    AbstractPolytope{N<:Real} <: AbstractPolyhedron{N}

Abstract type for compact convex polytopic sets.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

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

"""
    chebyshev_center(P::AbstractPolytope{N};
                     backend=default_polyhedra_backend(P, N),
                     compute_radius::Bool=false) where {N<:Real}

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
of a polytope.

### Input

- `P`              -- polytope
- `backend`        -- (optional, default: `default_polyhedra_backend(P, N)`) the
                      backend for polyhedral computations
- `compute_radius` -- (optional; default: `false`) option to additionally return
                      the radius of the largest ball enclosed by `P` around the
                      Chebyshev center

### Output

If `compute_radius` is `false`, the result is the Chebyshev center of `P`.
If `compute_radius` is `true`, the result is the pair `(c, r)` where `c` is the
Chebyshev center of `P` and `r` is the radius of the largest ball with center
`c` enclosed by `P`.

### Notes

The Chebyshev center is the center of a largest Euclidean ball enclosed by `P`.
In general, the center of such a ball is not unique (but the radius is).

Formally, letting ``∂P`` denote the boundary of ``P``, the Chebyshev center is

```math
    \\arg\\min_c \\min_{x ∈ ∂P} ‖x - c‖₂.
```
"""
function chebyshev_center(P::AbstractPolytope{N};
                          backend=default_polyhedra_backend(P, N),
                          compute_radius::Bool=false
                         ) where {N<:Real}
    require(:Polyhedra; fun_name="chebyshev_center")
    require(:JuMP; fun_name="chebyshev_center")
    require(:GLPK; fun_name="chebyshev_center")
    @assert applicable(polyhedron, P) "try converting the set to HPolytope " *
        "or VPolytope first"

    solver = JuMP.with_optimizer(GLPK.Optimizer)
    c, r = Polyhedra.chebyshevcenter(polyhedron(P; backend=backend), solver)

    if compute_radius
        return c, r
    end
    return c
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

# the "backend" argument is ignored, used for dispatch
function _vertices_list(P::AbstractPolytope, backend)
    return vertices_list(P)
end
