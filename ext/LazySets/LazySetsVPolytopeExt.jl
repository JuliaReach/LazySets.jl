using Base: extrema
using LazySets: LinearMapVRep, convex_hull!, dim, default_polyhedra_backend
using LazySets.HPolytopeModule: HPolytope
using LazySets.IntervalModule: Interval, _constraints_list_Vector
using LazySets.VPolygonModule: VPolygon
using LazySets.VPolytopeModule: VPolytope
using ReachabilityBase.Arrays: projection_matrix
using ReachabilityBase.Require: require
import LazySets.API: constraints_list, project, convex_hull
import LazySets: tohrep, _linear_map_vrep

"""
    constraints_list(P::VPolytope)

### Algorithm

For one- and two-dimensional sets, we respectively convert to an `Interval` or a
`VPolytope` and call the corresponding `constraints_list` function.
For higher-dimensional sets, we use `tohrep` to compute the constraint
representation and call the corresponding `constraints_list` function.
"""
function constraints_list(P::VPolytope)
    n = dim(P)
    if n == 1
        return _constraints_list_Vector(convert(Interval, P))
    elseif n == 2
        return constraints_list(convert(VPolygon, P))
    else
        return constraints_list(tohrep(P))
    end
end

# TODO move back to VPolytopeModule once `convex_hull!` is an API function
@inline function _linear_map_vrep(M::AbstractMatrix, P::VPolytope,
                                  algo::LinearMapVRep=LinearMapVRep(nothing);  # ignored
                                  apply_convex_hull::Bool=false)
    vlist = broadcast(v -> M * v, P.vertices)
    if apply_convex_hull
        convex_hull!(vlist)
    end
    return VPolytope(vlist)
end

@validate function project(P::VPolytope, block::AbstractVector{Int}; kwargs...)
    # the following cannot happen anymore because of `@validate`
    @assert !isempty(P.vertices) "empty VPolytope cannot be projected"

    m = length(block)
    if m == 1
        l, h = extrema(P, block[1])
        return Interval(l, h)
    end

    M = projection_matrix(block, dim(P), eltype(P.vertices))
    πvertices = broadcast(v -> M * v, P.vertices)

    if m == 2
        return VPolygon(πvertices; apply_convex_hull=true)
    else
        return VPolytope(πvertices)
    end
end

# TODO outsource to helper function in LazySets
"""
    convex_hull(P1::VPolytope, P2::VPolytope; [backend]=nothing)

Compute the convex hull of two polytopes in vertex representation.

### Input

- `P1`      -- polytope in vertex representation
- `P2`      -- polytope in vertex representation
- `backend` -- (optional, default: `nothing`) the polyhedral computation backend

### Output

The `VPolytope` obtained by the concrete convex hull of `P1` and `P2`.

### Notes

This function takes the union of the vertices of each polytope and then relies
on a concrete convex-hull algorithm.

The implementation relies on the polyhedral backend, which can be specified
using the `backend` keyword argument.

For performance reasons, it is suggested to use the `CDDLib.Library()` backend.
"""
@validate function convex_hull(P1::VPolytope, P2::VPolytope; backend=nothing)
    vunion = [P1.vertices; P2.vertices]
    convex_hull!(vunion; backend=backend)
    return VPolytope(vunion)
end

"""
    tohrep(P::VPolytope; [backend]=nothing)

Transform a polytope in vertex representation to a polytope in constraint
representation.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral computations; see
               [Polyhedra's documentation](https://juliapolyhedra.github.io/) for further
               information

### Output

A `HPolytope` as the constraint representation of `P`.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.

If `backend` is `nothing` and a backend is required, it defaults to `default_polyhedra_backend(P)`.
"""
function tohrep(P::VPolytope; backend=nothing)
    @assert !isempty(P.vertices) "cannot convert an empty polytope in vertex " *
                                 "representation to constraint representation"

    n = dim(P)
    if n == 1
        return HPolytope(_constraints_list_Vector(convert(Interval, P)))
    elseif n == 2
        return convert(HPolytope, convert(VPolygon, P))
    end

    require(LazySets, :Polyhedra; fun_name="tohrep")

    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    Q = LazySets.polyhedron(P; backend=backend)
    R = convert(HPolytope, Q)
    T = _output_type(P)  # help with type inference
    return R::T
end

function _output_type(::VPolytope)
    return HPolytope{N,Vector{N}} where {N}
end

function _output_type(::VPolytope{Float64})
    return HPolytope{Float64,Vector{Float64}}
end
