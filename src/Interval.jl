import IntervalArithmetic
import IntervalArithmetic: AbstractInterval
import Base:+, -, *, ∈, ⊆

export Interval,
       dim, σ, center,
       low, high, vertices_list

"""
    Interval{N<:Real, IN <: AbstractInterval{N}} <: AbstractPointSymmetricPolytope{N}

Type representing an interval on the real line. Mathematically, it is of the
form

```math
[a, b] := \\{ a ≤ x ≤ b \\} ⊆ \\mathbb{R}.
```

### Fields

- `dat` -- data container for the given interval

### Notes

This type relies on the [IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/)
library for representation of intervals and arithmetic operations.

### Examples

Unidimensional intervals are symbolic representations of a real closed interval.
This type requires the user to load the `IntervalArithmetic` library, since
artithmetic operations rely on that module.

We can create intervals in different ways, the simpler way is to pass a pair
of numbers:

```jldoctest interval_constructor
julia> x = Interval(0.0, 1.0)
LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```
or a 2-vector:

```jldoctest interval_constructor
julia> x = Interval([0.0, 1.0])
LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```

Note that if the package `IntervalArithmetic` is loaded in the current scope,
you have to prepend the `LazySets` to the interval type, since there is
a name conflict otherwise.

```jldoctest interval_constructor
julia> using IntervalArithmetic
WARNING: using IntervalArithmetic.Interval in module Main conflicts with an existing identifier.

julia> x = LazySets.Interval(IntervalArithmetic.Interval(0.0, 1.0))
LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])

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
LazySets.Interval{Rational{Int64},IntervalArithmetic.AbstractInterval{Rational{Int64}}}([0//1, 2//1])
```
"""
struct Interval{N<:Real, IN <: AbstractInterval{N}} <: AbstractPointSymmetricPolytope{N}
   dat::IN
end
# type-less convenience constructor
Interval(interval::IN) where {N, IN <: AbstractInterval{N}} = Interval{N, IN}(interval)

# constructor that takes two numbers
Interval(lo::N, hi::N) where {N} = Interval(IntervalArithmetic.Interval(lo, hi))

Interval(lo::Rational{N}, hi::Rational{N}) where {N} = Interval{Rational{N}, IntervalArithmetic.AbstractInterval{Rational{N}}}(IntervalArithmetic.Interval(lo, hi))

# constructor from a vector
Interval(x::AbstractVector{N}) where {N} = Interval(IntervalArithmetic.Interval(x[1], x[2]))

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
    σ(d::V, x::Interval{N, IN})::V where {N, IN <: AbstractInterval{N}, V<:AbstractVector{N}}

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Support vector in the given direction.
"""
function σ(d::V, x::Interval{N, IN})::V where {N, IN <: AbstractInterval{N}, V<:AbstractVector{N}}
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
    center(x::Interval)

Return the interval's center.

### Input

- `x` -- interval

### Output

The center, or midpoint, of ``x``.
"""
⊆(x::Interval, y::Interval) = x.dat ⊆ y.dat

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
    vertices_list(x::Interval)

Return the list of vertices of this interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval represented as two one-dimensional vectors.
"""
vertices_list(x::Interval) = [[low(x)], [high(x)]]
