"""
    reflect(x::AbstractVector, H::Hyperplane)

Reflect (mirror) a vector in a hyperplane.

### Input

- `x` -- point/vector
- `H` -- hyperplane

### Output

The reflection of `x` in `H`.

### Algorithm

The reflection of a point ``x`` in the hyperplane ``a ⋅ x = b`` is

```math
    x − 2 \\frac{x ⋅ a − b}{a ⋅ a} a
```

where ``u · v`` denotes the dot product.
"""
@commutative function reflect(x::AbstractVector, H::Hyperplane)
    return _reflect_point_hyperplane(x, H.a, H.b)
end

function _reflect_point_hyperplane(x, a, b)
    return x - 2 * (dot(x, a) - b) / dot(a, a) * a
end
