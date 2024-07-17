function project(H::Hyperplane{N}, block::AbstractVector{Int}; kwargs...) where {N}
    if constrained_dimensions(H) ⊆ block
        return Hyperplane(H.a[block], H.b)
    else
        return Universe{N}(length(block))
    end
end

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
