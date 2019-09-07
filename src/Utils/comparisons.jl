# Some functions in this file are inspired from Polyhedra.jl
#
export set_rtol, set_ztol, set_atol

"""
    Tolerance{N<:Number}

Type that represents the tolerances for a given numeric type.

### Fields

- `rtol` -- relative tolerance
- `ztol` -- zero tolerance or absolute tolerance for comparison against zero
- `atol` -- absolute tolerance

### Notes

The type `Tolerance`, parametric in the numeric type `N`, is used to store default
values for numeric comparisons. It is mutable and setting the value of a field
affects the getter functions hence it can be used to fix the tolerance globally
in `LazySets`.

Default values are defined for the most commonly used numeric types, and for those
cases when other numeric types are needed one can extend the default values
as explained next.

The cases `Float64` and `Rational` are special in the sense that they are the most
commonly used types in applications. Getting and setting default tolerances
is achieved with the functions `_rtol` and `set_rtol` (and similarly for the other
tolerances); the implementation creates an instance of `Tolerance{Float64}`
(resp. `Tolerance{Rational}`) and sets some default values. Again since `Tolerance`
is mutable, setting a value is possible e.g. `set_rtol(Type{Float64}, ε)` for some
floating-point `ε`.

For all other cases, a dictionary mapping numeric types to instances of `Tolerance`
for that numeric type is used. For floating-point types, a default value has been
defined through `default_tolerance` as follows:

```julia
default_tolerance(N::Type{<:AbstractFloat}) = Tolerance(Base.rtoldefault(N), sqrt(eps(N)), zero(N))
```
Hence to set a single tolerance (either `rtol`, `ztol` or `atol`) for a given
floating-point type, use the corresponding `set_rtol` function, while the values
which have not been set will be pulled from `default_tolerance`. If you would like
to define the three default values at once, or are computing with a non floating-point
numeric type, you can just extend `default_tolerance(N::Type{<:Number})`.
"""
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

# global default tolerances for other numeric types
TOL_N = Dict{Type{<:Number}, Tolerance}()
_rtol(N::Type{<:Number}) = get!(TOL_N, N, default_tolerance(N)).rtol
_ztol(N::Type{<:Number}) = get!(TOL_N, N, default_tolerance(N)).ztol
_atol(N::Type{<:Number}) = get!(TOL_N, N, default_tolerance(N)).atol

set_rtol(N::Type{NT}, ε::NT) where {NT<:Number} = begin
    if N ∉ keys(TOL_N)
        TOL_N[N] = default_tolerance(N)
    end
    TOL_N[N].rtol = ε
end

set_ztol(N::Type{NT}, ε::NT) where {NT<:Number} = begin
    if N ∉ keys(TOL_N)
        TOL_N[N] = default_tolerance(N)
    end
    TOL_N[N].ztol = ε
end

set_atol(N::Type{NT}, ε::NT) where {NT<:Number} = begin
    if N ∉ keys(TOL_N)
        TOL_N[N] = default_tolerance(N)
    end
    TOL_N[N].atol = ε
end

default_tolerance(N::Type{<:Number}) = error("default tolerance for numeric type $N is not defined")
default_tolerance(N::Type{<:AbstractFloat}) = Tolerance(Base.rtoldefault(N), sqrt(eps(N)), zero(N))

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

"""
    minmax(A, B, C)

Compute the minimum and maximum of three numbers A, B, C.

### Input

- `A` -- first number
- `B` -- second number
- `C` -- third number

### Output

The minimum and maximum of the three inputed numbers.

### Examples

```jldoctest
julia> LazySets.minmax(1.4, 52.4, -5.2)
(-5.2, 52.4)
```
"""
function minmax(A, B, C)
    if A > B
        min, max = B, A
    else
        min, max = A, B
    end
    if C > max
        max = C
    elseif C < min
        min = C
    end
    return min, max
end

"""
    arg_minmax(A, B, C)

Compute the index (1, 2, 3) of the minimum and maximum of three numbers A, B, C.

### Input

- `A` -- first number
- `B` -- second number
- `C` -- third number

### Output

The index of the minimum and maximum of the three inputed numbers.

### Examples

```jldoctest
julia> LazySets.arg_minmax(1.4, 52.4, -5.2)
(3, 2)
```
"""
function arg_minmax(A, B, C)
    if A > B
        min, max = B, A
        imin, imax = 2, 1
    else
        min, max = A, B
        imin, imax = 1, 2
    end
    if C > max
        imax = 3
    elseif C < min
        imin = 3
    end
    return imin, imax
end
