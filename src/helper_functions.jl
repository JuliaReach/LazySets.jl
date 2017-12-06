import Base.<=

export unit_step,
       jump2pi

"""
    unit_step(x::Real)

The unit step function, which returns 1 iff x is greater than or equal to zero.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``0`` otherwise.

### Notes

This function can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> unit_step.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
unit_step(x::Real) = ifelse(x < zero(x), oftype(one(x),-1), one(x))

"""
    jump2pi(x::Float64)

Return ``x + 2π`` if ``x`` is negative, otherwise return ``x``.

### Input

- `x` -- real scalar

### Output

``x + 2π`` if ``x`` is negative, ``x`` otherwise.

### Examples

```jldoctest
julia> jump2pi(0.0)
0.0
julia> jump2pi(-0.5)
5.783185307179586
julia> jump2pi(0.5)
0.5
```
"""
jump2pi(x::Float64) = x < 0 ? 2.0 * pi + x : x

"""
    u <= v

Compares two 2D vectors by their direction.

### Input

- `u` --  first 2D direction
- `v` --  second 2D direction

### Output

True iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).
"""
function <=(u::AbstractVector{Float64}, v::AbstractVector{Float64})
    return jump2pi(atan2(u[2], u[1])) <= jump2pi(atan2(v[2], v[1]))
end
