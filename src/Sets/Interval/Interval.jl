"""
    Interval{N} <: AbstractHyperrectangle{N}

Type representing an interval on the real line.
Mathematically, it is of the form

```math
[a, b] := \\{ a ≤ x ≤ b \\} ⊆ ℝ.
```

### Fields

- `dat` -- data container for the given interval

### Notes

This type relies on the
[IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/)
library for representation of intervals and arithmetic operations.

### Examples

Unidimensional intervals are symbolic representations of a real closed interval.

We can create intervals in different ways. The simplest way is to pass a pair
of numbers:

```jldoctest interval_constructor
julia> x = Interval(0.0, 1.0)
Interval{Float64}([0, 1])
```
A 2-vector is also possible:

```jldoctest interval_constructor
julia> x = Interval([0.0, 1.0])
Interval{Float64}([0, 1])
```

An interval can also be constructed from an `IntervalArithmetic.Interval`.
Note that if the package `IntervalArithmetic` is loaded in the current scope,
you have to prepend the package names to `Interval`s since there is a name
conflict otherwise.

```jldoctest interval_constructor
julia> using IntervalArithmetic

julia> x = LazySets.Interval(IntervalArithmetic.interval(0.0, 1.0))
LazySets.IntervalModule.Interval{Float64}([0, 1])

julia> dim(x)
1

julia> center(x)
1-element Vector{Float64}:
 0.5
```

The usual pairwise arithmetic operators `-` and `*` are interpreted in the
standard way known in interval arithmetic, so the results are `Interval`s.
Note that `+` is generally used for the lazy Minkowksi sum in this library.

Intervals of other numeric types can be created as well, e.g., a rational
interval:

```jldoctest
julia> Interval(0//1, 2//1)
Interval{Rational{Int64}}([0//1, 2//1])
```
"""
struct Interval{N} <: AbstractHyperrectangle{N}
    dat::IA.Interval{N}

    function Interval(dat::IA.Interval{N}) where {N}
        @assert isfinite(IA.inf(dat)) && isfinite(IA.sup(dat)) "intervals must be bounded"

        return new{N}(dat)
    end
end

# constructor from two numbers with type promotion
function Interval(lo::N1, hi::N2) where {N1,N2}
    N = promote_type(N1, N2)
    return Interval(IA.interval(N(lo), N(hi)))
end

# constructor from a vector
function Interval(x::AbstractVector)
    @assert length(x) == 2 "cannot construct an Interval from a " *
                           "$(length(x))-dimensional vector"
    return @inbounds Interval(x[1], x[2])
end

# constructor from a single number
function Interval(x::Real)
    return Interval(IA.interval(x))
end
