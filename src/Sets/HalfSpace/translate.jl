"""
# Extended help

    translate(hs::HalfSpace, v::AbstractVector; [share]::Bool=false)

### Notes

The normal vector of the half-space (vector `a` in `a⋅x ≤ b`) is shared with
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
