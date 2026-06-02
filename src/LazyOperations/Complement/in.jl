"""
    in(x::AbstractVector, C::Complement)

Check whether a given point is contained in the complement of a set.

### Input

- `x` -- point/vector
- `C` -- complement of a set

### Output

`true` iff the vector is contained in the complement.

### Algorithm

```math
    x ∈ X^C ⟺ x ∉ X
```
"""
@validate function in(x::AbstractVector, C::Complement)
    return x ∉ C.X
end
