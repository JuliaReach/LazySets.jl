"""
    ∈(x::AbstractVector, hs::HalfSpace)

Check whether a given point is contained in a half-space.

### Input

- `x` -- point/vector
- `hs` -- half-space

### Output

`true` iff ``x ∈ hs``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x ≤ b``.
"""
function ∈(x::AbstractVector, hs::HalfSpace)
    return _leq(dot(x, hs.a), hs.b)
end
