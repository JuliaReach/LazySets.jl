import IntervalArithmetic
using IntervalArithmetic: AbstractInterval
import Base: +, -, *, ∈, ⊆, rand, min, max

export Interval,
       dim, σ, center,
       vertices_list,
       isflat

"""
    Interval{N<:Real, IN<:AbstractInterval{N}} <: AbstractHyperrectangle{N}

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
you have to prepend `LazySets.` to the `Interval` type since there is a name
conflict otherwise.

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
Interval{Rational{Int64},AbstractInterval{Rational{Int64}}}([0//1, 2//1])
```
"""
struct Interval{N<:Real, IN<:AbstractInterval{N}} <: AbstractHyperrectangle{N}
    dat::IN
end

# convenience constructor without type parameter for Rational
Interval(interval::IN) where {N<:Rational, IN<:AbstractInterval{N}} =
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

Return the support vector of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{N}, x::Interval{N}) where {N<:Real}
    @assert length(d) == dim(x)
    return d[1] > zero(N) ? high(x) : low(x)
end

"""
    ρ(d::AbstractVector{N}, x::Interval{N}) where {N<:Real}

Evaluate the support function of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector{N}, x::Interval{N}) where {N<:Real}
    @assert length(d) == dim(x)
    return d[1] * (d[1] > zero(N) ? max(x) : min(x))
end

"""
    center(x::Interval{N})::Vector{N} where {N<:Real}

Return the interval's center.

### Input

- `x` -- interval

### Output

The center, or midpoint, of ``x``.
"""
function center(x::Interval{N})::Vector{N} where {N<:Real}
    return [IntervalArithmetic.mid(x.dat)]
end

"""
    +(x::Interval{N}, y::Interval{N}) where {N<:Real}

Return the sum of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The sum of the intervals as a new `Interval` set.
"""
function +(x::Interval{N}, y::Interval{N}) where {N<:Real}
    return Interval(x.dat + y.dat)
end

"""
    -(x::Interval{N}, y::Interval{N}) where {N<:Real}

Return the difference of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The difference of the intervals as a new `Interval` set.
"""
function -(x::Interval{N}, y::Interval{N}) where {N<:Real}
    return Interval(x.dat - y.dat)
end

"""
```
    *(x::Interval{N}, y::Interval{N}) where {N<:Real}
```

Return the product of the intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

The product of the intervals as a new `Interval` set.
"""
function *(x::Interval{N}, y::Interval{N}) where {N<:Real}
    return Interval(x.dat * y.dat)
end

"""
    ∈(v::AbstractVector{N}, x::Interval{N}) where {N<:Real})

Return whether a vector is contained in the interval.

### Input

- `v` -- one-dimensional vector
- `x` -- interval

### Output

`true` iff `x` contains `v`'s first component.
"""
function ∈(v::AbstractVector{N}, x::Interval{N}) where {N<:Real}
    return v[1] ∈ x.dat
end

"""
    ∈(v::N, x::Interval{N}) where {N<:Real}

Return whether a number is contained in the interval.

### Input

- `v` -- scalar
- `x` -- interval

### Output

`true` iff `x` contains `v`.
"""
function ∈(v::N, x::Interval{N}) where {N<:Real}
    return v ∈ x.dat
end

"""
    min(x::Interval{N})::N where {N<:Real}

Return the lower component of an interval.

### Input

- `x` -- interval

### Output

The lower (`lo`) component of the interval.
"""
function min(x::Interval{N})::N where {N<:Real}
    return x.dat.lo
end

"""
    max(x::Interval{N})::N where {N<:Real}

Return the higher or upper component of an interval.

### Input

- `x` -- interval

### Output

The higher (`hi`) component of the interval.
"""
function max(x::Interval{N})::N where {N<:Real}
    return x.dat.hi
end

"""
    low(x::Interval{N})::Vector{N} where {N<:Real}

Return the lower coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the lower coordinate of the interval.
"""
function low(x::Interval{N})::Vector{N} where {N<:Real}
    return N[x.dat.lo]
end

"""
    high(x::Interval{N})::Vector{N} where {N<:Real}

Return the higher coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the higher coordinate of the interval.
"""
function high(x::Interval{N})::Vector{N} where {N<:Real}
    return N[x.dat.hi]
end

"""
    an_element(x::Interval{N})::Vector{N} where {N<:Real}

Return some element of an interval.

### Input

- `x` -- interval

### Output

The left border (`min(x)`) of the interval.
"""
function an_element(x::Interval{N})::Vector{N} where {N<:Real}
    return [min(x)]
end

"""
    rand(::Type{Interval}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Interval{N}

Create a random interval.

### Input

- `Interval` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 1) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

A random interval.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Interval};
              N::Type{<:Real}=Float64,
              dim::Int=1,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Interval{N}
    @assert dim == 1 "cannot create a random Interval of dimension $dim"
    rng = reseed(rng, seed)
    x = randn(rng, N)
    y = randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

"""
    vertices_list(x::Interval{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of this interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval represented as two one-dimensional vectors.
"""
function vertices_list(x::Interval{N})::Vector{Vector{N}} where {N<:Real}
    return [[min(x)], [max(x)]]
end

"""
    translate(x::Interval{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) an interval by a given vector.

### Input

- `x` -- interval
- `v` -- translation vector

### Output

A translated interval.

### Algorithm

We add the vector to the left and right of the interval.
"""
function translate(x::Interval{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(x) "cannot translate a $(dim(x))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Interval(x.dat + v[1])
end


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(x::Interval{N}, i::Int)::N where {N<:Real}

Return the box radius of an interval in a given dimension.

### Input

- `x` -- interval
- `i` -- dimension index (must be `1`)

### Output

The box radius in the given dimension.
"""
function radius_hyperrectangle(x::Interval{N}, i::Int)::N where {N<:Real}
    @assert i == 1 "an interval is one-dimensional"
    return (max(x) - min(x)) / N(2)
end

"""
    radius_hyperrectangle(x::Interval{N})::Vector{N} where {N<:Real}

Return the box radius of an interval in every dimension.

### Input

- `x` -- interval

### Output

The box radius of the interval (a one-dimensional vector).
"""
function radius_hyperrectangle(x::Interval{N})::Vector{N} where {N<:Real}
    return [radius_hyperrectangle(x, 1)]
end

"""
    isflat(I::Interval)::Bool

Determine whether an interval is flat, i.e. whether its extreme values coincide.

### Input

- `I` -- interval

### Output

A boolean which is `true` if the interval is flat and `false` otherwise.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` when applied to the diameter of the interval.
Hence, this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(I::Interval)::Bool
    return isapproxzero(IntervalArithmetic.diam(I.dat))
end

"""
    plot_recipe(I::Interval{N}, [ε]::N=zero(N)) where {N<:Real}

Convert an interval to a pair `(x, y)` of points for plotting.

### Input

- `I` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.

### Notes

We consider the interval as a line segment with y coordinate equal to zero.
"""
function plot_recipe(I::Interval{N}, ε::N=zero(N)) where {N<:Real}
    return [min(I), max(I)], zeros(N, 2)
end
