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

    require(@__MODULE__, :LazySets; fun_name="tohrep")

    n = dim(P)
    if n == 1
        Q = convert(Interval, P)
        return convert(HPolytope, Q)
    elseif n == 2
        Q = convert(VPolygon, P)
        return convert(HPolytope, Q)
    end

    require(@__MODULE__, :Polyhedra; fun_name="tohrep")

    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    return convert(HPolytope, LazySets.polyhedron(P; backend=backend))
end
