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
function translate(P::HPoly, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end
