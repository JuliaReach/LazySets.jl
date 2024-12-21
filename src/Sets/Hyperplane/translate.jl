"""
# Extended help

    translate(H::Hyperplane, v::AbstractVector; share::Bool=false)

### Notes

The normal vector of the hyperplane (vector ``a`` in ``a⋅x = b``) is shared with
the original hyperplane if `share == true`.

### Algorithm

A hyperplane ``a⋅x = b`` is transformed to the hyperplane ``a⋅x = b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(H::Hyperplane, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(H) "cannot translate a $(dim(H))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    a = share ? H.a : copy(H.a)
    b = H.b + dot(H.a, v)
    return Hyperplane(a, b)
end
