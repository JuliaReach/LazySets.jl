"""
    sign_cadlag(x::Real)

This function works like the sign function but is ``1`` for input ``0``.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``-1`` otherwise.

### Notes

This is the sign function right-continuous at zero (see
[càdlàg function](https://en.wikipedia.org/wiki/C%C3%A0dl%C3%A0g)).
It can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> LazySets.sign_cadlag.([-0.6, 1.3, 0.0])
3-element Vector{Float64}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::Real)
    return x < zero(x) ? -one(x) : one(x)
end
