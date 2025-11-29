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
    require(@__MODULE__, :LazySets; fun_name="tovrep")

    n = dim(P)
    if n == 1
        if isempty(P)
            N = eltype(P)
            return VPolytope(Vector{N}[])
        else
            Q = convert(Interval, P)
            return convert(VPolytope, Q)
        end
    elseif n == 2
        Q = convert(HPolygon, P)
        return convert(VPolytope, Q)
    end

    require(@__MODULE__, :Polyhedra; fun_name="tovrep")

    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    return convert(VPolytope, LazySets.polyhedron(P; backend=backend))
end
