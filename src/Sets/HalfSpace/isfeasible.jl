
"""
    isfeasible(constraints::AbstractVector{<:HalfSpace}, [witness]::Bool=false;
               [solver]=nothing)

Check for feasibility of a list of linear constraints.

### Input

- `constraints` -- list of linear constraints
- `witness`     -- (optional; default: `false`) flag for witness production
- `solver`      -- (optional; default: `nothing`) LP solver

### Output

If `witness` is `false`, the result is a `Bool`.

If `witness` is `true`, the result is a pair `(res, w)` where `res` is a `Bool`
and `w` is a witness point/vector.

### Algorithm

This implementation converts the constraints to matrix-vector form via
`tosimplehrep` and then calls `isfeasible` on the result.
"""
function isfeasible(constraints::AbstractVector{<:HalfSpace},
                    witness::Bool=false; solver=nothing)
    # the caller should verify that there is at least one constraint
    A, b = tosimplehrep(constraints)
    return isfeasible(A, b, witness; solver=solver)
end
