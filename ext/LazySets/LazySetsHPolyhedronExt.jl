import LazySets

using LazySets.EmptySetModule: EmptySet
using LazySets.HPolygonModule: HPolygon
using LazySets.HPolyhedronModule: HPoly
using LazySets.IntervalModule: Interval
using LazySets.VPolytopeModule: VPolytope
using LazySets: default_polyhedra_backend, dim, remove_redundant_constraints!
using ReachabilityBase.Require: require
import LazySets: remove_redundant_constraints, tovrep

"""
    tovrep(P::HPoly; [backend]=nothing)

Transform a polytope in constraint representation to a polytope in vertex
representation.

### Input

- `P`       -- polytope in constraint representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral computations

### Output

A `VPolytope` which is a vertex representation of the given polytope in
constraint representation.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
For further information on the supported backends see [Polyhedra's
documentation](https://juliapolyhedra.github.io/).

If `backend` is `nothing` and a backend is required, it defaults to `default_polyhedra_backend(P)`.
"""
function tovrep(P::HPoly; backend=nothing)
    n = dim(P)
    if n == 1
        if isempty(P)
            N = eltype(P)
            return VPolytope(Vector{N}[])
        end
        Q = convert(Interval, P)
    elseif n == 2
        Q = convert(HPolygon, P)
    else
        require(LazySets, :Polyhedra; fun_name="tovrep")

        if isnothing(backend)
            backend = default_polyhedra_backend(P)
        end

        Q = LazySets.polyhedron(P; backend=backend)
    end
    return convert(VPolytope, Q)
end

"""
    remove_redundant_constraints(P::HPoly{N}; [backend]=nothing) where {N}

Remove the redundant constraints in a polyhedron in constraint representation.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Output

A polyhedron equivalent to `P` but with no redundant constraints, or an empty
set if `P` is detected to be empty, which may happen if the constraints are
infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:HalfSpace})` for details.
"""
function remove_redundant_constraints(P::HPoly; backend=nothing)
    Pred = copy(P)
    if remove_redundant_constraints!(Pred; backend=backend)
        return Pred
    else # the polyhedron P is empty
        N = eltype(P)
        return EmptySet{N}(dim(P))
    end
end
