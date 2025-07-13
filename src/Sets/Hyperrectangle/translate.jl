"""
# Extended help

    translate(H::Hyperrectangle, v::AbstractVector; [share]::Bool=false)

### Input

- `H`     -- hyperrectangle
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Notes

The radius vector is shared with the original hyperrectangle if `share == true`.
"""
@validate function translate(H::Hyperrectangle, v::AbstractVector; share::Bool=false)
    c = center(H) + v
    radius = share ? H.radius : copy(H.radius)
    return Hyperrectangle(c, radius)
end

function translate!(H::Hyperrectangle, v::AbstractVector)
    @assert length(v) == dim(H) "cannot translate a $(dim(H))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    H.center .+= v
    return H
end
