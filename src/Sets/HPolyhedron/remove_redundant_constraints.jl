# see ext/LazySets/LazySetsHPolyhedronExt.jl

"""
    remove_redundant_constraints!(P::HPoly{N}; [backend]=nothing) where {N}

Remove the redundant constraints of a polyhedron in constraint representation;
the polyhedron is updated in-place.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Output

`true` if the method was successful and the polyhedron `P` is modified by
removing its redundant constraints, and `false` if `P` is detected to be empty,
which may happen if the constraints are infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:HalfSpace})` for details.
"""
function remove_redundant_constraints!(P::HPoly; backend=nothing)
    return remove_redundant_constraints!(P.constraints; backend=backend)
end
