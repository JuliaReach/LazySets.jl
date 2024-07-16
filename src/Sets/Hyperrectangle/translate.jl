"""
    translate(H::Hyperrectangle, v::AbstractVector; [share]::Bool=false)

Translate (i.e., shift) a hyperrectangle by a given vector.

### Input

- `H`     -- hyperrectangle
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated hyperrectangle.

### Notes

The radius vector is shared with the original hyperrectangle if `share == true`.

### Algorithm

We add the vector to the center of the hyperrectangle.
"""
function translate(H::Hyperrectangle, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(H) "cannot translate a $(dim(H))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
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
