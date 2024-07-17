"""
    ∈(x::AbstractVector, H::Hyperplane)

Check whether a given point is contained in a hyperplane.

### Input

- `x` -- point/vector
- `H` -- hyperplane

### Output

`true` iff ``x ∈ H``.

### Algorithm

We just check whether ``x`` satisfies ``a⋅x = b``.
"""
function ∈(x::AbstractVector, H::Hyperplane)
    return _isapprox(dot(H.a, x), H.b)
end
