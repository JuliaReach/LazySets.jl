"""
    translate(L::Line2D, v::AbstractVector; [share]::Bool=false)

Translate (i.e., shift) a 2D line by a given vector.

### Input

- `L`     -- 2D line
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated line.

### Notes

The normal vector of the line (vector `a` in `a⋅x = b`) is shared with the
original line if `share == true`.

### Algorithm

A line ``a⋅x = b`` is transformed to the line ``a⋅x = b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(L::Line2D, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    a = share ? L.a : copy(L.a)
    b = L.b + dot(L.a, v)
    return Line2D(a, b)
end
