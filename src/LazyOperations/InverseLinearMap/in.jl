"""
    in(x::AbstractVector, ilm::InverseLinearMap)

Check whether a given point is contained in the inverse linear map of a set.

### Input

- `x`   -- point/vector
- `ilm` -- inverse linear map of a set

### Output

`true` iff ``x ∈ ilm``.

### Algorithm

This implementation does not explicitly invert the matrix since it uses the
property ``x ∈ M^{-1}⋅X`` iff ``M⋅x ∈ X``.

### Examples

```jldoctest
julia> ilm = LinearMap([0.5 0.0; 0.0 -0.5], BallInf([0., 0.], 1.));

julia> [1.0, 1.0] ∈ ilm
false

julia> [0.1, 0.1] ∈ ilm
true
```
"""
@validate function in(x::AbstractVector, ilm::InverseLinearMap)
    y = ilm.M * x
    return y ∈ ilm.X
end
