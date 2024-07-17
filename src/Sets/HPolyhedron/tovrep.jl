"""
    tovrep(P::HPoly; [backend]=default_polyhedra_backend(P))

Transform a polytope in constraint representation to a polytope in vertex
representation.

### Input

- `P`       -- polytope in constraint representation
- `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
               backend for polyhedral computations

### Output

A `VPolytope` which is a vertex representation of the given polytope in
constraint representation.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
For further information on the supported backends see [Polyhedra's
documentation](https://juliapolyhedra.github.io/).
"""
function tovrep(P::HPoly; backend=default_polyhedra_backend(P))
    require(@__MODULE__, :LazySets; fun_name="tovrep")
    require(@__MODULE__, :Polyhedra; fun_name="tovrep")

    P = polyhedron(P; backend=backend)
    return VPolytope(P)
end
