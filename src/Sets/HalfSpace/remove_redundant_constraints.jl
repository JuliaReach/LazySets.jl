"""
    remove_redundant_constraints(constraints::AbstractVector{<:HalfSpace}; backend=nothing)

Remove the redundant constraints of a given list of linear constraints.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `nothing`) the backend used to solve the
                   linear program
### Output

The list of constraints with the redundant ones removed, or an empty list if the
constraints are infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `[remove_redundant_constraints!(::AbstractVector{<:HalfSpace})](@ref)` for
details.
"""
function remove_redundant_constraints(constraints::AbstractVector{<:HalfSpace}; backend=nothing)
    constraints_copy = copy(constraints)
    if remove_redundant_constraints!(constraints_copy; backend=backend)
        return constraints_copy
    else  # the constraints are infeasible
        return Vector{eltype(constraints)}(undef, 0)
    end
end

"""
    remove_redundant_constraints!(constraints::AbstractVector{<:HalfSpace}; [backend]=nothing)

Remove the redundant constraints of a given list of linear constraints; the list
is updated in-place.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `nothing`) the backend used to solve the
                   linear program
### Output

`true` if the removal was successful and the list of constraints `constraints`
is modified by removing the redundant constraints, and `false` only if the
constraints are infeasible.

### Notes

Note that the result may be `true` even if the constraints are infeasible.
For example, ``x ≤ 0 && x ≥ 1`` will return `true` without removing any
constraint.
To check if the constraints are infeasible, use
`isempty(HPolyhedron(constraints))`.

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

If there are `m` constraints in `n` dimensions, this function checks one by one
if each of the `m` constraints is implied by the remaining ones.

To check if the `k`-th constraint is redundant, an LP is formulated using the
constraints that have not yet been removed.
If, at an intermediate step, it is detected that a subgroup of the constraints
is infeasible, this function returns `false`.
If the calculation finished successfully, this function returns `true`.

For details, see [Fukuda's Polyhedra
FAQ](https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node24.html).
"""
function remove_redundant_constraints!(constraints::AbstractVector{<:HalfSpace};
                                       backend=nothing)
    if isempty(constraints)
        return true
    end
    N = eltype(first(constraints))
    if isnothing(backend)
        backend = default_lp_solver(N)
    end
    A, b = tosimplehrep(constraints)
    m, n = size(A)
    non_redundant_indices = 1:m

    i = 1  # counter over reduced (= non-redundant) constraints

    for j in 1:m  # loop over original constraints
        α = A[j, :]
        Ar = A[non_redundant_indices, :]
        br = b[non_redundant_indices]
        @assert br[i] == b[j]
        br[i] += one(N)
        lp = linprog(-α, Ar, '<', br, -Inf, Inf, backend)
        if is_lp_infeasible(lp.status)
            # the polyhedron is empty
            return false
        elseif is_lp_optimal(lp.status)
            objval = -lp.objval
            if _leq(objval, b[j])
                # the constraint is redundant
                non_redundant_indices = setdiff(non_redundant_indices, j)
            else
                # the constraint is not redundant
                i += 1
            end
        else
            throw(ArgumentError("LP is not optimal; the status of the LP is $(lp.status)"))
        end
    end

    deleteat!(constraints, setdiff(1:m, non_redundant_indices))
    return true
end
