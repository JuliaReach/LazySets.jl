"""
# Extended help

    translate(P::HPoly, v::AbstractVector; [share]::Bool=false)

### Input

- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.
"""
@validate function translate(P::HPoly, v::AbstractVector; share::Bool=false)
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

"""
# Extended help

    translate!(P::HPoly, v::AbstractVector)

### Algorithm

We translate every constraint in place.
"""
@validate function translate!(P::HPoly, v::AbstractVector)
    constraints = P.constraints
    @inbounds for i in eachindex(constraints)
        constraints[i] = translate(constraints[i], v)
    end
    return P
end
