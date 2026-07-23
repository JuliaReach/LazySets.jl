"""
    remove_redundant_vertices(P::VPolytope; [backend]=nothing, [solver]=nothing)

Return the polytope obtained by removing the redundant vertices of the given
polytope in vertex representation.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral
               computations; see `default_polyhedra_backend(P)` or
               [Polyhedra's documentation](https://juliapolyhedra.github.io/)
               for further information
- `solver`  -- (optional, default: `nothing`) the linear programming
               solver used in the backend, if needed; see
               `default_lp_solver_polyhedra(N)`

### Output

A new polytope such that its vertices are the convex hull of the given polytope.

### Notes

The optimization problem associated to removing redundant vertices is handled
by `Polyhedra`.
If the polyhedral computations backend requires an LP solver, but it has not
been specified, we use `default_lp_solver_polyhedra(N)` to define such solver.
Otherwise, the redundancy-removal function of the polyhedral backend is used.
"""
function remove_redundant_vertices(P::VPolytope; backend=nothing, solver=nothing)
    return _remove_redundant_vertices(P; backend, solver)
end

# see ext/PolyhedraExt.jl
function _remove_redundant_vertices(P; backend=nothing, solver=nothing)
    mod = Base.get_extension(@__MODULE__, :PolyhedraExt)
    require(mod, :Polyhedra; fun_name="remove_redundant_vertices")
    error()
end
