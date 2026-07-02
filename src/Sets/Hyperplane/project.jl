# see ext/LazySets/LazySetsHyperplaneExt.jl

"""
    project(x::AbstractVector, H::Hyperplane)

Project a point onto a hyperplane.

### Input

- `x` -- point
- `H` -- hyperplane

### Output

The projection of `x` onto `H`.

### Algorithm

The projection of ``x`` onto the hyperplane of the form ``a⋅x = b`` is

```math
    x - \\dfrac{a (a⋅x - b)}{‖a‖²}
```
"""
function project(x::AbstractVector, H::Hyperplane)
    return x - H.a * (dot(H.a, x) - H.b) / norm(H.a, 2)^2
end
