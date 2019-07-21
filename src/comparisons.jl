# Some functions in this file are inspired from Polyhedra.jl
#
export set_rtol

# struct to contain the tolerances for a given numberic type
mutable struct Tolerance{N<:Number}
    rtol::N
    ztol::N
    atol::N
end

# global Float64 tolerances
const _TOL_F64 = Tolerance(Base.rtoldefault(Float64), Float64(10)*sqrt(eps(Float64)), zero(Float64))

_rtol(N::Type{Float64}) = _TOL_F64.rtol
_ztol(N::Type{Float64}) = _TOL_F64.ztol
_atol(N::Type{Float64}) = _TOL_F64.atol

set_rtol(N::Type{Float64}, ε::Float64) = _TOL_F64.rtol = ε
set_ztol(N::Type{Float64}, ε::Float64) = _TOL_F64.ztol = ε
set_atol(N::Type{Float64}, ε::Float64) = _TOL_F64.atol = ε

# global rational tolerances
const _TOL_RAT = Tolerance(zero(Rational), zero(Rational), zero(Rational))

_rtol(N::Type{<:Rational}) = _TOL_RAT.rtol
_ztol(N::Type{<:Rational}) = _TOL_RAT.ztol
_atol(N::Type{<:Rational}) = _TOL_RAT.atol

set_rtol(N::Type{<:Rational}, ε::Rational) = _TOL_RAT.rtol = ε
set_ztol(N::Type{<:Rational}, ε::Rational) = _TOL_RAT.ztol = ε
set_atol(N::Type{<:Rational}, ε::Rational) = _TOL_RAT.atol = ε

# global default tolerances (cannot be set)
_rtol(N::Type{<:AbstractFloat}) = Base.rtoldefault(N)
_ztol(N::Type{<:AbstractFloat}) = sqrt(eps(N))
_atol(N::Type{<:AbstractFloat}) = zero(N)

"""
    _leq(x::N, y::N; [kwargs...]) where {N<:Real}

Determine if `x` is smaller than or equal to `y`.

### Input

- `x`      -- number
- `y`      -- another number (of the same numeric type as `x`)
- `kwargs` -- not used

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

This is a fallback implementation for numbers of type `Real`. If the arguments
are floating point numbers, see `_leq(x::AbstractFloat, y::AbstractFloat)`.
"""
_leq(x::N, y::N; kwargs...) where {N<:Real} = x <= y

"""
    _leq(x::N, y::M; [kwargs...]) where {N<:Real, M<:Real}

Determine if `x` is smaller than or equal to `y`.

### Input

- `x`      -- number
- `y`      -- another number (of possibly different numeric type than `x`)
- `kwargs` -- optional arguments; see `?_leq` for the available options 

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

This implementation calls Julia's `promote(x, y)` function, which converts all
arguments to a common numeric type, returning them as a tuple. The conversion
is such that the common type to which the values are converted can represent
them as faithfully as possible.
"""
_leq(x::N, y::M; kwargs...) where {N<:Real, M<:Real} =
    _leq(promote(x, y)...; kwargs...)

"""
    _geq(x::Real, y::Real; [kwargs...])

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
_geq(x::Real, y::Real; kwargs...) = _leq(y, x; kwargs...)

"""
    isapproxzero(x::N; ztol::Real=_ztol(N)) where {N<:Real}

Determine if `x` is approximately zero.

### Input

- `x`    -- number
- `ztol` -- (optional, default: `_ztol(N)`) tolerance against zero

### Output

A boolean that is `true` iff `x ≈ 0`.

### Algorithm

It is considered that `x ≈ 0` whenever `x` (in absolute value) is smaller than
the tolerance for zero, `ztol`.
"""
function isapproxzero(x::N; ztol::Real=_ztol(N)) where {N<:Real}
    return abs(x) <= ztol
end

"""
    _isapprox(x::N, y::N;
              rtol::Real=_rtol(N),
              ztol::Real=_ztol(N),
              atol::Real=_atol(N)) where {N<:Real}

Determine if `x` is approximately equal to `y`.

### Input

- `x`    -- number
- `y`    -- another number (of the same numeric type as `x`)
- `rtol` -- (optional, default: `_rtol(N)`) relative tolerance
- `ztol` -- (optional, default: `_ztol(N)`) absolute tolerance for comparison
            against zero
- `atol` -- (optional, default: `_atol(N)`) absolute tolerance

### Output

A boolean that is `true` iff `x ≈ y`.

### Algorithm

We first check if `x` and `y` are both approximately zero, using
`isapproxzero(x, y)`.
If that fails, we check if `x ≈ y`, using Julia's `isapprox(x, y)`.
In the latter check we use `atol` absolute tolerance and `rtol` relative
tolerance.

Comparing to zero with default tolerances is a special case in Julia's
`isapprox`, see the last paragraph in `?isapprox`. This function tries to
combine `isapprox` with its default values and a branch for `x ≈ y ≈ 0` which
includes `x == y == 0` but also admits a tolerance `ztol`.

Note that if `x = ztol` and `y = -ztol`, then `|x-y| = 2*ztol` and still
`_isapprox` returns `true`.
"""
function _isapprox(x::N, y::N;
                   rtol::Real=_rtol(N),
                   ztol::Real=_ztol(N),
                   atol::Real=_atol(N)) where {N<:Real}
    if isapproxzero(x, ztol=ztol) && isapproxzero(y, ztol=ztol)
        return true
    else
        return isapprox(x, y, rtol=rtol, atol=atol)
    end
end

# generic "dense"
function _isapprox(x::AbstractVector{N}, y::AbstractVector{N};
                   rtol::Real=_rtol(N),
                   ztol::Real=_ztol(N),
                   atol::Real=_atol(N)) where {N<:Real}
    n = length(x)
    if length(x) != length(y)
        return false
    end
    @inbounds for i in 1:n
        if !_isapprox(x[i], y[i], rtol=rtol, ztol=ztol, atol=atol)
            return false
        end
    end
    return true
end

# sparse
function _isapprox(x::SparseVector{N}, y::SparseVector{N};
                   rtol::Real=_rtol(N),
                   ztol::Real=_ztol(N),
                   atol::Real=_atol(N)) where {N<:Real}
    @assert length(x) == length(y)
    return x.nzind == y.nzind && _isapprox(x.nzval, y.nzval, rtol=rtol, ztol=ztol, atol=atol)
end

"""
    ispermutation(u::AbstractVector{T}, v::AbstractVector)::Bool where {T}

Check that two vectors contain the same elements up to reordering.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`true` iff the vectors are identical up to reordering.

### Examples

```jldoctest
julia> LazySets.ispermutation([1, 2, 2], [2, 2, 1])
true

julia> LazySets.ispermutation([1, 2, 2], [1, 1, 2])
false
```

### Notes

Containment check is performed using `LazySets._in(e, v)`, so in the case of
floating point numbers, the precision to which the check is made is determined
by the type of elements in `v`. See `_in` and `_isapprox` for more information.
"""
function ispermutation(u::AbstractVector{T}, v::AbstractVector)::Bool where {T}
    if length(u) != length(v)
        return false
    end
    occurrence_map = Dict{T, Int}()
    has_duplicates = false
    for e in u
        if !_in(e, v)
            return false
        end
        if haskey(occurrence_map, e)
            occurrence_map[e] += 1
            has_duplicates = true
        else
            occurrence_map[e] = 1
        end
    end
    if has_duplicates
        for e in v
            if !haskey(occurrence_map, e) || occurrence_map[e] == 0
                return false
            end
            occurrence_map[e] -= 1
        end
    end
    return true
end

# alias for Julia's containment check 
function _in(x, itr) where {T}
    return x ∈ itr
end

# approximate containment check for numbers in floating point 
function _in(x::AbstractVector{T}, itr) where {T<:AbstractFloat}
    return any(y -> _isapprox(x, y), itr)
end

"""
    _leq(x::N, y::N;
         rtol::Real=_rtol(N),
         ztol::Real=_ztol(N),
         atol::Real=_atol(N)) where {N<:AbstractFloat}

Determine if `x` is smaller than or equal to `y`.

### Input

- `x`    -- number
- `y`    -- another number (of the same numeric type as `x`)
- `rtol` -- (optional, default: `_rtol(N)`) relative tolerance
- `ztol` -- (optional, default: `_ztol(N)`) absolute tolerance for comparison
            against zero
- `atol` -- (optional, default: `_atol(N)`) absolute tolerance

### Output

A boolean that is `true` iff `x <= y`.

### Algorithm

The `x <= y` comparison is split into `x < y` or `x ≈ y`; the latter is
implemented by extending Juila's built-in `isapprox(x, y)` with an absolute
tolerance that is used to compare against zero.
"""
function _leq(x::N, y::N;
              rtol::Real=_rtol(N),
              ztol::Real=_ztol(N),
              atol::Real=_atol(N)) where {N<:AbstractFloat}
    return x <= y || _isapprox(x, y, rtol=rtol, ztol=ztol, atol=atol)
end
