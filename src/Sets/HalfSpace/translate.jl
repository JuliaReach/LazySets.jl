"""
    translate(hs::HalfSpace, v::AbstractVector; [share]::Bool=false)

Translate (i.e., shift) a half-space by a given vector.

### Input

- `hs`    -- half-space
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated half-space.

### Notes

The normal vectors of the half-space (vector `a` in `a⋅x ≤ b`) is shared with
the original half-space if `share == true`.

### Algorithm

A half-space ``a⋅x ≤ b`` is transformed to the half-space ``a⋅x ≤ b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(hs::HalfSpace, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(hs) "cannot translate a $(dim(hs))-dimensional " *
                                 "set by a $(length(v))-dimensional vector"
    a = share ? hs.a : copy(hs.a)
    b = hs.b + dot(hs.a, v)
    return HalfSpace(a, b)
end
