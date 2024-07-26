"""
    vertices_list(P::HPolytope; [backend]=nothing, [prune]::Bool=true)

Return a list of the vertices of a polytope in constraint representation.

### Input

- `P`       -- polytope in constraint representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral
               computations
- `prune`   -- (optional, default: `true`) flag to remove redundant vertices

### Output

A list of the vertices.

### Algorithm

If the polytope is one-dimensional (resp. two-dimensional), it is converted to
an interval (resp. polygon in constraint representation) and then the respective
optimized `vertices_list` implementation is used.

It is possible to use the `Polyhedra` backend in the one- and two-dimensional
case as well by passing a `backend`.

If the polytope is not two-dimensional, the concrete polyhedra-manipulation
library `Polyhedra` is used. The actual computation is performed by a given
backend; for the default backend used in `LazySets` see
`default_polyhedra_backend(P)`. For further information on the supported
backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function vertices_list(P::HPolytope; backend=nothing, prune::Bool=true)
    N = eltype(P)
    if isempty(P.constraints)
        return Vector{N}(Vector{N}(undef, 0))  # illegal polytope
    elseif isnothing(backend)
        # use efficient implementations in 1D and 2D
        n = dim(P)
        if n == 1
            return vertices_list_1d(P)
        elseif n == 2
            require(@__MODULE__, :LazySets; fun_name="vertices_list")

            return vertices_list(convert(HPolygon, P; prune=prune))
        end
    end

    require(@__MODULE__, :Polyhedra; fun_name="vertices_list")
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end
    Q = Polyhedra.polyhedron(P; backend=backend)
    if prune
        Polyhedra.removevredundancy!(Q; ztol=_ztol(N))
    end
    return collect(Polyhedra.points(Q))
end
