"""
# Extended help

    translate(v::AbstractVector, P::HPolygon; [share]::Bool=false)

### Input

- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.

### Algorithm

We translate every constraint.
"""
function translate(P::HPolygon, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    return HPolygon(constraints; sort_constraints=false,
                    check_boundedness=false, prune=false)
end
