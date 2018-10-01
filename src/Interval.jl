import IntervalArithmetic
import IntervalArithmetic: AbstractInterval
import Base:+, -, *, ∈, ⊆

export Interval,
       dim, σ, center,
       low, high, vertices_list

"""
    Interval{N<:Real, IN <: AbstractInterval{N}} <: AbstractHyperrectangle{N}

Type representing an interval on the real line.
Mathematically, it is of the form

```math
[a, b] := \\{ a ≤ x ≤ b \\} ⊆ \\mathbb{R}.
```

### Fields

- `dat` -- data container for the given interval

### Notes

This type relies on the
[IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/)
library for representation of intervals and arithmetic operations.

### Examples

Unidimensional intervals are symbolic representations of a real closed interval.

We can create intervals in different ways, the simpler way is to pass a pair
of numbers:

```jldoctest interval_constructor
julia> x = Interval(0.0, 1.0)
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```
or a 2-vector:

```jldoctest interval_constructor
julia> x = Interval([0.0, 1.0])
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```

Note that if the package `IntervalArithmetic` is loaded in the current scope,
you have to prepend the `LazySets` to the interval type, since there is
a name conflict otherwise.

```jldoctest interval_constructor
julia> using IntervalArithmetic
WARNING: using IntervalArithmetic.Interval in module Main conflicts with an existing identifier.

julia> x = Interval(IntervalArithmetic.Interval(0.0, 1.0))
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])

julia> dim(x)
1

julia> center(x)
1-element Array{Float64,1}:
 0.5
```

This type is such that the usual pairwise arithmetic operators `+`, `-`, `*` trigger
the corresponding interval arithmetic backend method, and return a new
`Interval` object. For the symbolic Minkowksi sum, use `MinkowskiSum` or `⊕`.

Interval of other numeric types can be created as well, eg. a rational interval:

```jldoctest interval_constructor
julia> Interval(0//1, 2//1)
Interval{Rational{Int64},IntervalArithmetic.AbstractInterval{Rational{Int64}}}([0//1, 2//1])
```
"""
struct Interval{N<:Real, IN <: AbstractInterval{N}} <: AbstractHyperrectangle{N}
   dat::IN
end

@static if VERSION < v"0.7-"
    # convenience constructor without type parameter
    Interval(interval::IN) where {N<:Real, IN <: AbstractInterval{N}} =
        Interval{N, IN}(interval)
end

# convenience constructor without type parameter for Rational
Interval(interval::IN) where {N<:Rational, IN <: AbstractInterval{N}} =
    Interval{N, IntervalArithmetic.AbstractInterval{N}}(interval)

# constructor from two numbers
Interval(lo::N, hi::N) where {N<:Real} =
    Interval(IntervalArithmetic.Interval(lo, hi))

# constructor from two rational numbers
Interval(lo::N, hi::N) where {N<:Rational} =
    Interval{N, IntervalArithmetic.AbstractInterval{N}}(
        IntervalArithmetic.Interval(lo, hi))

# constructor from a vector
function Interval(x::AbstractVector{N}) where {N<:Real}
    @assert length(x) == 2 "vector for Interval constructor has to be 2D"
    Interval(IntervalArithmetic.Interval(x[1], x[2]))
end

# constructor from a rational vector
function Interval(x::AbstractVector{N}) where {N<:Rational}
    @assert length(x) == 2 "vector for Interval constructor has to be 2D"
    Interval{N, IntervalArithmetic.AbstractInterval{N}}(
        IntervalArithmetic.Interval(x[1], x[2]))
end

"""
    dim(x::Interval)::Int

Return the ambient dimension of an interval.

### Input

- `x` -- interval

### Output

The integer 1.
"""
dim(x::Interval)::Int = 1

"""
    σ(d::AbstractVector{N}, x::Interval{N}) where {N<:Real}

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{N}, x::Interval{N}) where {N<:Real}
    @assert length(d) == dim(x)
    return d[1] > 0 ? [x.dat.hi] : [x.dat.lo]
end

"""
    center(x::Interval)

Return the interval's center.

### Input

- `x` -- interval

### Output

The center, or midpoint, of ``x``.
"""
center(x::Interval) = [IntervalArithmetic.mid(x.dat)]

"""
    +(x::Interval, y::Interval)

Return the sum of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The sum of the intervals as a new `Interval` set.
"""
+(x::Interval, y::Interval) = Interval(x.dat + y.dat)

"""
    -(x::Interval, y::Interval)

Return the difference of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The difference of the intervals as a new `Interval` set.
"""
-(x::Interval, y::Interval) = Interval(x.dat - y.dat)

"""
```
    *(x::Interval, y::Interval)
```

Return the product of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The product of the intervals as a new `Interval` set.
"""
*(x::Interval, y::Interval) = Interval(x.dat * y.dat)

"""
    ∈(v::AbstractVector, x::Interval)

Return whether a vector is contained in the interval.

### Input

- `v` -- one-dimensional vector
- `x` -- interval

### Output

`true` iff `x` contains `v`'s first component.
"""
∈(v::AbstractVector, x::Interval) = v[1] ∈ x.dat

"""
    ∈(v::N, x::Interval) where {N}

Return whether a number is contained in the interval.

### Input

- `v` -- scalar
- `x` -- interval

### Output

`true` iff `x` contains `v`.
"""
∈(v::N, x::Interval) where {N} = v ∈ x.dat

"""
    low(x::Interval)

Return the lower component of an interval.

### Input

- `x` -- interval

### Output

The lower (`lo`) component of the interval.
"""
low(x::Interval) = x.dat.lo

"""
    high(x::Interval)

Return the higher or upper component of an interval.

### Input

- `x` -- interval

### Output

The higher (`hi`) component of the interval.
"""
high(x::Interval) = x.dat.hi

"""
    an_element(x::Interval{N})::Vector{N} where {N<:Real}

Return some element of an interval.

### Input

- `x` -- interval

### Output

The left border (`low(x)`) of the interval.
"""
function an_element(x::Interval{N})::Vector{N} where {N<:Real}
    return [low(x)]
end

"""
    vertices_list(x::Interval)

Return the list of vertices of this interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval represented as two one-dimensional vectors.
"""
vertices_list(x::Interval) = [[low(x)], [high(x)]]


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(I::Interval{N}, i::Int)::N where {N<:Real}

Return the box radius of an interval in a given dimension.

### Input

- `I` -- interval
- `i` -- dimension index (must be `1`)

### Output

The box radius in the given dimension.
"""
function radius_hyperrectangle(I::Interval{N}, i::Int)::N where {N<:Real}
    @assert i == 1 "an interval is one-dimensional"
    return high(x) - low(x)
end

"""
    radius_hyperrectangle(I::Interval{N})::Vector{N} where {N<:Real}

Return the box radius of an interval in every dimension.

### Input

- `I` -- interval

### Output

The box radius of the interval (a one-dimensional vector).
"""
function radius_hyperrectangle(I::Interval{N})::Vector{N} where {N<:Real}
    return [high(x) - low(x)]
end
