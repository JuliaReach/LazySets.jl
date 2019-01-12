const ABSZTOL(N::Type{<:AbstractFloat}) = sqrt(eps(N))
const ABSZTOL(N::Type{Rational{INNER}}) where {INNER} = zero(N)

"""
    _leq(x::N, y::N) where {N<:Real}

Return `true` if the real `x` is smaller or equal than `y` and false otherwise.

### Input

- `x` -- number
- `y` -- another number (of the same numeric type as `x`)

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

This is a fallback implementation for numbers of type `Real`. If the arguments
are floating point numbers, see `_leq(x::AbstractFloat, y::AbstractFloat)`.
"""
_leq(x::N, y::N) where {N<:Real} = x <= y

"""
    _leq(x::N, y::M) where {N<:Real, M<:Real}

Determine if `x` is smaller than or equal to `y`.

### Input

- `x` -- number
- `y` -- another number (of possibly different numeric type than `x`)

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

This implementation calls Julia's `promote(x, y)` function, which converts all
arguments to a common numeric type, returning them as a tuple. The conversion
is such that the common type to which the values are converted can represent
them as faithfully as possible.
"""
_leq(x::N, y::M) where {N<:Real, M<:Real} = _leq(promote(x, y)...)

"""
    _geq(x::N, y::N) where {N<:Real}

Determine if `x` is greater than or equal to `y`.

### Input

- `x` -- number
- `y` -- another number (of the same numeric type as `x`)

### Output

A boolean that is `true` iff `x >= y`.

### Algorithm

This function falls back to `_leq(y, x)`. See the documentation of `_leq` for
further details. 
"""
_geq(x::N, y::N) where {N<:Real} = _leq(y, x)

"""
    _geq(x::N, y::M) where {N<:Real, M<:Real}

Determine if `x` is greater than or equal to `y`.

### Input

- `x` -- number
- `y` -- another number (of possibly different numeric type than `x`)

### Output

A boolean that is `true` iff `x >= y`.

### Algorithm

This function falls back to `_leq(y, x)`, with type promotion if needed. See the
documentation of `_leq` for further details.
"""
_geq(x::N, y::M) where {N<:Real, M<:Real} = _geq(promote(x, y)...)

"""
    isapproxzero(x::N; kwargs...) where {N<:Real}

Determine if `x` is approximately zero.

### Input

- `x`      -- number
- `kwargs` -- ignored
 
### Output

A boolean that is `true` iff `x ≈ 0`.

### Algorithm

This is a fallback implementation for any real `x` such that `x ≈ 0` is `true`
whenever `x` is equal to zero in the same numeric type as `x`.
"""
function isapproxzero(x::N; kwargs...) where {N<:Real}
    return x == zero(N)
end

"""
    isapproxzero(x::N; ztol::Real=ABSZTOL(N)) where {N<:AbstractFloat}

Returns `true` if `x` is approximately zero `false` otherwise.

### Input

- `x`    -- number
- `ztol` -- (optional, default: `ABSZTOL`) tolerance against zero

### Output

A boolean that is `true` iff `x ≈ 0`.

### Algorithm

It is considered that `x ≈ 0` whenever `x` (in absolute value) is smaller than the
tolerance for zero, `ztol`.
"""
function isapproxzero(x::N; ztol::Real=ABSZTOL(N)) where {N<:AbstractFloat}
    return abs(x) < ztol
end

"""
    _isapprox(x::N, y::N;
              rtol::Real=Base.rtoldefault(N),
              ztol::Real=ABSZTOL(N),
              atol::Real=zero(N)) where {N<:AbstractFloat}

Determine if `x` is approximately equal to `y`.

### Input

- `x`    -- number
- `y`    -- another number (of the same numeric type than `x`)
- `rtol` -- (optional, default: `Base.rtoldefault(N)`) relative tolerance
- `ztol` -- (optional, default: `ABSZTOL(N)`) absolute tolerance for comparison against zero
- `atol` -- absolute tolerance

### Output

A boolean that is `true` iff `x ≈ y`.

### Algorithm

This comparison is split into both `x` and `y` approximately zero checked using
`isapproxzero(x, y)`, or `x ≈ y` checked using Julia's `isapprox(x, y)`. In the
latter function we use zero absolute tolerance and `rtol` relative tolerance.
"""
function _isapprox(x::N, y::N;
                   rtol::Real=Base.rtoldefault(N),
                   ztol::Real=ABSZTOL(N),
                   atol::Real=zero(N)) where {N<:AbstractFloat}
    if isapproxzero(x, ztol=ztol) && isapproxzero(y, ztol=ztol)
        return true
    else
        return isapprox(x, y, rtol=rtol, atol=zero(N))
    end
end

"""
    _leq(x::N, y::N;
         rtol::Real=Base.rtoldefault(N),
         ztol::Real=ABSZTOL(N),
         atol::Real=zero(N)) where {N<:AbstractFloat}

Returns `true` if the real `x` is smaller or equal than `y` and false otherwise.

### Input

- `x`    -- number
- `y`    -- another number (of the same numeric type than `x`)
- `rtol` -- (optional, default: `Base.rtoldefault(N)`) relative tolerance
- `ztol` -- (optional, default: `ABSZTOL(N)`) absolute tolerance for comparison against zero
- `atol` -- absolute tolerance

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

The `x <= y` comparison is split into `x < y` or `x ≈ y`; the latter is implemented
by extending Juila's built-in `isapprox(x, y)` with an absolute tolerance that is
used to compare against zero.
"""
function _leq(x::N, y::N;
              rtol::Real=Base.rtoldefault(N),
              ztol::Real=ABSZTOL(N),
              atol::Real=zero(N)) where {N<:AbstractFloat}
    return x < y || _isapprox(x, y, rtol=rtol, ztol=ztol, atol=atol)
end
