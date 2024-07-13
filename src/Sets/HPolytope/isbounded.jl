"""
    isbounded(P::HPolytope, [use_type_assumption]::Bool=true)

Determine whether a polytope in constraint representation is bounded.

### Input

- `P`                   -- polytope in constraint representation
- `use_type_assumption` -- (optional, default: `true`) flag for ignoring the
                           type assumption that polytopes are bounded

### Output

`true` if `use_type_assumption` is activated.
Otherwise, `true` iff `P` is bounded.
"""
function isbounded(P::HPolytope, use_type_assumption::Bool=true)
    if use_type_assumption
        return true
    end
    return isbounded(P.constraints)
end
