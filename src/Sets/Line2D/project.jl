# see ext/LazySets/LazySetsLine2DExt.jl

"""
    project(x::AbstractVector, L::Line2D)

Project a point onto a 2D line.

### Input

- `x` -- point/vector
- `L` -- 2D line

### Output

The projection of `x` onto `L`.

### Algorithm

The projection of ``x`` onto a line of the form ``a⋅x = b`` is

```math
    x - \\dfrac{a (a⋅x - b)}{‖a‖²}.
```
"""
function project(x::AbstractVector, L::Line2D)
    return x - L.a * (dot(L.a, x) - L.b) / norm(L.a, 2)^2
end
