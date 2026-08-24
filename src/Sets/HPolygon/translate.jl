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
@validate function translate(P::HPolygon, v::AbstractVector; share::Bool=false)
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    return HPolygon(constraints; sort_constraints=false,
                    check_boundedness=false, prune=false)
end

"""
# Extended help

    translate!(P::HPolygon, v::AbstractVector)

### Algorithm

We translate every constraint in place. Since translation preserves the
normal vectors (only the offsets change), the constraints remain sorted.
"""
@validate function translate!(P::HPolygon, v::AbstractVector)
    constraints = P.constraints
    @inbounds for i in eachindex(constraints)
        constraints[i] = translate(constraints[i], v)
    end
    return P
end
