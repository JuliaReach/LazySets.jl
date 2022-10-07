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

"""
    minmax(a, b, c)

Compute the minimum and maximum of three numbers a, b, c.

### Input

- `a` -- number
- `b` -- number
- `c` -- number

### Output

The minimum and maximum of three given numbers.

### Examples

```jldoctest
julia> LazySets.minmax(1.4, 52.4, -5.2)
(-5.2, 52.4)
```
"""
function minmax(a, b, c)
    if a > b
        min, max = b, a
    else
        min, max = a, b
    end
    if c > max
        max = c
    elseif c < min
        min = c
    end
    return min, max
end

"""
    arg_minmax(a, b, c)

Compute the indices of the minimum and maximum of three numbers a, b, c.

### Input

- `a` -- first number
- `b` -- second number
- `c` -- third number

### Output

The indices of the minimum and maximum of the three given numbers.

### Examples

```jldoctest
julia> LazySets.arg_minmax(1.4, 52.4, -5.2)
(3, 2)
```
"""
function arg_minmax(a, b, c)
    if a > b
        min, max = b, a
        imin, imax = 2, 1
    else
        min, max = a, b
        imin, imax = 1, 2
    end
    if c > max
        imax = 3
    elseif c < min
        imin = 3
    end
    return imin, imax
end
