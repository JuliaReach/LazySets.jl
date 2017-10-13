"""
    unit_step(x)

The unit step function, which returns 1 if and only if x is greater or equal
than zero.

### Input

- ``x`` -- a floating point number

NOTE:

This function can be used with vector-valued arguments via the dot operator.

EXAMPLES:

    julia> unit_step([-0.6, 1.3, 0.0])
    3-element Array{Float64,1}:
     -1.0
     1.0
     1.0
"""
unit_step(x::Float64) = ifelse(x < 0, oftype(one(x),-1), one(x))

"""
    jump2pi(x)

Return `x + 2*pi` and only if `x` is negative.

### Input

- ``x`` -- a floating point number

EXAMPLES:

    julia> jump2pi(0.0)
    0.0
    julia> jump2pi(-0.5)
    5.783185307179586
    julia> jump2pi(0.5)
    0.5
"""
jump2pi(x::Float64) = x < 0 ? 2.0 * pi + x : x

export unit_step, jump2pi