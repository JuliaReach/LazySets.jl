import Base: +, -, *, ∈, ⊆, rand, min, max

export Interval,
       dim, σ, center,
       vertices_list,
       isflat,
       linear_map
constraints_list

"""
    Interval{N} <: AbstractHyperrectangle{N}

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
WARNING: using IntervalArithmetic.Interval in module Main conflicts with an existing identifier.

julia> x = LazySets.Interval(IntervalArithmetic.Interval(0.0, 1.0))
Interval{Float64}([0, 1])

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

```jldoctest interval_constructor
julia> Interval(0//1, 2//1)
Interval{Rational{Int64}}([0//1, 2//1])
```
"""
struct Interval{N} <: AbstractHyperrectangle{N}
    dat::IA.Interval{N}

    function Interval(dat::IA.Interval{N}) where {N}
        @assert isfinite(dat.lo) && isfinite(dat.hi) "intervals must be bounded"

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

isoperationtype(::Type{<:Interval}) = false

"""
    dim(x::Interval)

Return the ambient dimension of an interval.

### Input

- `x` -- interval

### Output

``1``.
"""
dim(x::Interval) = 1

"""
    σ(d::AbstractVector, x::Interval)

Return the support vector of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, x::Interval)
    @assert length(d) == dim(x) "a $(length(d))-dimensional vector is " *
                                "incompatible with an $(dim(x))-dimensional set"
    N = promote_type(eltype(d), eltype(x))
    return @inbounds d[1] > zero(N) ? high(x) : low(x)
end

"""
    ρ(d::AbstractVector, x::Interval)

Evaluate the support function of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, x::Interval)
    @assert length(d) == dim(x) "a $(length(d))-dimensional vector is " *
                                "incompatible with an $(dim(x))-dimensional set"
    N = promote_type(eltype(d), eltype(x))
    return @inbounds d[1] * (d[1] > zero(N) ? max(x) : min(x))
end

"""
    center(x::Interval)

Return the center of an interval.

### Input

- `x` -- interval

### Output

The center, or midpoint, of `x`.
"""
function center(x::Interval)
    return [_center(x)]
end

"""
    center(H::Interval, i::Int)

Return the center along a given dimension of a interval.

### Input

- `x` -- interval
- `i` -- dimension of interest

### Output

The center along a given dimension of the interval.
"""
@inline function center(x::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return _center(x)
end

# returns a number, not a vector
function _center(x::Interval)
    return IA.mid(x.dat)
end

"""
    -(x::Interval, y::Interval)

Return the difference of two intervals (in the interval-arithmetic sense).

### Input

- `x` -- interval
- `y` -- interval

### Output

The difference of the intervals as a new `Interval` set.
"""
function -(x::Interval, y::Interval)
    return Interval(x.dat - y.dat)
end

"""
```
    *(x::Interval, y::Interval)
```

Return the product of two intervals (in the interval-arithmetic sense).

### Input

- `x` -- interval
- `y` -- interval

### Output

The product of the intervals as a new `Interval` set.
"""
function *(x::Interval, y::Interval)
    return Interval(x.dat * y.dat)
end

"""
    ∈(v::AbstractVector, x::Interval))

Check whether a given point is contained in an interval.

### Input

- `v` -- point/vector
- `x` -- interval

### Output

`true` iff ``v ∈ x``.
"""
function ∈(v::AbstractVector, x::Interval)
    @assert length(v) == 1 "a $(length(v))-dimensional vector is " *
                           "incompatible with an interval"
    return @inbounds v[1] ∈ x.dat
end

"""
    ∈(v::Number, x::Interval)

Check whether a number is contained in an interval.

### Input

- `v` -- scalar
- `x` -- interval

### Output

`true` iff `x` contains the singleton `[v]`.
"""
function ∈(v::Number, x::Interval)
    return v ∈ x.dat
end

"""
    min(x::Interval)

Return the lower component of an interval.

### Input

- `x` -- interval

### Output

The lower (`lo`) component of the interval (a number).
"""
function min(x::Interval)
    return x.dat.lo
end

"""
    max(x::Interval)

Return the higher or upper component of an interval.

### Input

- `x` -- interval

### Output

The higher (`hi`) component of the interval (a number).
"""
function max(x::Interval)
    return x.dat.hi
end

"""
    low(x::Interval)

Return the lower coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the lower coordinate of the interval.
"""
function low(x::Interval)
    return [x.dat.lo]
end

"""
    high(x::Interval)

Return the higher coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the higher coordinate of the interval.
"""
function high(x::Interval)
    return [x.dat.hi]
end

"""
    an_element(x::Interval)

Return some element of an interval.

### Input

- `x` -- interval

### Output

The left border (`low(x)`) of the interval.
"""
function an_element(x::Interval)
    return low(x)
end

"""
    rand(::Type{Interval}; [N]::Type{<:Real}=Float64, [dim]::Int=1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 1 "cannot create a random Interval of dimension $dim"
    rng = reseed!(rng, seed)
    x = randn(rng, N)
    y = randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

"""
    vertices_list(x::Interval; kwargs...)

Return the list of vertices of an interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval, which are two one-dimensional vectors,
or just one if the interval is degenerate (the endpoints match within the
working tolerance).
"""
function vertices_list(x::Interval; kwargs...)
    a = min(x)
    b = max(x)
    return _isapprox(a, b) ? [[a]] : [[a], [b]]
end

"""
    constraints_list(x::Interval)

Return the list of constraints of an interval.

### Input

- `x` -- interval

### Output

The list of constraints of the interval represented as two one-dimensional
half-spaces.
"""
function constraints_list(x::Interval)
    N = eltype(x)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2)
    e₁ = SingleEntryVector(1, 1, one(N))
    @inbounds constraints[1] = HalfSpace(e₁, max(x))
    @inbounds constraints[2] = HalfSpace(-e₁, -min(x))
    return constraints
end

"""
    translate(x::Interval, v::AbstractVector)

Translate (i.e., shift) an interval by a given vector.

### Input

- `x` -- interval
- `v` -- translation vector

### Output

A translated interval.

### Algorithm

We add the vector to the left and right of the interval.

### Notes

An in-place version is not available because the `IntervalArithmetic.Interval`
type is immutable.
"""
function translate(x::Interval, v::AbstractVector)
    @assert length(v) == dim(x) "cannot translate a $(dim(x))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Interval(x.dat + @inbounds v[1])
end

function _radius(x::Interval{N}) where {N}
    return (max(x) - min(x)) / N(2)
end

"""
    radius_hyperrectangle(x::Interval{N}, i::Int) where {N}

Return the box radius of an interval in a given dimension.

### Input

- `x` -- interval
- `i` -- dimension index (must be `1`)

### Output

The box radius in the given dimension.
"""
function radius_hyperrectangle(x::Interval{N}, i::Int) where {N}
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return _radius(x)
end

"""
    radius_hyperrectangle(x::Interval)

Return the box radius of an interval in every dimension.

### Input

- `x` -- interval

### Output

The box radius of the interval (a one-dimensional vector).
"""
function radius_hyperrectangle(x::Interval)
    return [_radius(x)]
end

"""
    isflat(x::Interval)

Check whether an interval is flat, i.e., whether its extreme values coincide.

### Input

- `x` -- interval

### Output

`true` if the interval is flat and `false` otherwise.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` when applied to the diameter of the interval.
In other words, this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(x::Interval)
    return isapproxzero(IA.diam(x.dat))
end

"""
    plot_recipe(x::Interval{N}, [ε]=zero(N)) where {N}

Convert an interval to a pair `(x, y)` of points for plotting.

### Input

- `x` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.

### Notes

We consider the interval as a line segment with y coordinate equal to zero.
"""
function plot_recipe(x::Interval{N}, ε=zero(N)) where {N}
    return [min(x), max(x)], zeros(N, 2)
end

"""
    linear_map(M::AbstractMatrix, x::Interval)

Concrete linear map of an interval.

### Input

- `M` -- matrix
- `x` -- interval

### Output

Either an interval or a zonotope, depending on the leading dimension (i.e., the
number of rows) of `M`:

- If `size(M, 1) == 1`, the output is an `Interval` obtained by scaling `x` by
  the matrix `M`.
- If `size(M, 1) ≠ 1`, the output is a `Zonotope` with center `M * center(x)`
  and the single generator `M * g`, where `g = (high(x)-low(x))/2`.
"""
function linear_map(M::AbstractMatrix, x::Interval)
    @assert size(M, 2) == 1 "a linear map of size $(size(M)) " *
                            "cannot be applied to an interval"
    nout = size(M, 1)
    if nout == 1
        return _linear_map_interval(M, x)
    else
        return _linear_map_zonotope(M, x)
    end
end

function _linear_map_interval(M::AbstractMatrix, x::Interval)
    α = @inbounds M[1, 1]
    return Interval(α * x.dat)
end

function _linear_map_zonotope(M::AbstractMatrix, x::Interval)
    nout = size(M, 1)
    cx = _center(x)
    gx = cx - min(x)
    N = promote_type(eltype(M), eltype(x))
    c = Vector{N}(undef, nout)
    gen = Matrix{N}(undef, nout, 1)
    @inbounds for i in 1:nout
        c[i] = M[i, 1] * cx
        gen[i] = M[i, 1] * gx
    end
    return Zonotope(c, gen)
end

"""
    scale(α::Real, x::Interval)

Concrete scaling of an interval.

### Input

- `α` -- scalar
- `x` -- interval

### Output

The interval obtained by scaling the given interval.
"""
function scale(α::Real, x::Interval)
    return Interval(α * x.dat)
end

"""
    rectify(x::Interval{N}) where {N}

Concrete rectification of an interval.

### Input

- `x` -- interval

### Output

The `Interval` that corresponds to the rectification of `x`.

### Notes

Note that the result is an `Interval` even if the set becomes a singleton (which
is the case if the original interval was nonpositive).
"""
function rectify(x::Interval{N}) where {N}
    if x.dat.lo >= zero(N)
        # interval is already nonnegative
        return x
    else
        # lower end is negative
        return Interval(zero(N), max(x.dat.hi, zero(N)))
    end
end

"""
    diameter(x::Interval, [p]::Real=Inf)

Compute the diameter of an interval, defined as ``\\Vert b - a\\Vert`` in the
``p`-norm, where ``a`` (resp. ``b``) are the minimum (resp. maximum) value of
the interval.

### Input

- `x` -- interval
- `p` -- (optional, default: `Inf`) norm (ignored)

### Output

A real number representing the diameter.

### Notes

In one dimension all p-norms are identical.
"""
function diameter(x::Interval, p::Real=Inf)
    return IA.diam(x.dat)
end

"""
    split(x::Interval, k)

Partition an interval into `k` uniform sub-intervals.

### Input

- `x` -- interval
- `k` -- number of sub-intervals, possibly wrapped in a vector of length 1

### Output

A list of `k` `Interval`s.
"""
function split(x::Interval, k::AbstractVector{Int})
    @assert length(k) == 1 "an interval can only be split along one dimension"
    return split(x, @inbounds k[1])
end

function split(x::Interval, k::Int)
    @assert k > 0 "can only split into a positive number of intervals"
    return [Interval(x2) for x2 in mince(x.dat, k)]
end

"""
    ngens(x::Interval)

Return the number of generators of an interval.

### Input

- `x` -- interval

### Output

The number of generators.

### Algorithm

An interval has either one generator, or zero generators if it is a degenerated
interval of diameter zero.
"""
function ngens(x::Interval)
    return _isapprox(min(x), max(x)) ? 0 : 1
end

"""
    chebyshev_center_radius(x::Interval; [kwargs]...)

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of an interval.

### Input

- `x`      -- interval
- `kwargs` -- further keyword arguments (ignored)

### Output

The pair `(c, r)` where `c` is the Chebyshev center of `x` and `r` is the radius
of the largest ball with center `c` enclosed by `x`.

### Notes

The Chebyshev center of an interval is just the center of the interval.
"""
function chebyshev_center_radius(x::Interval; kwargs...)
    return center(x), _radius(x)
end

"""
    reflect(x::Interval)

Concrete reflection of a interval `x`, resulting in the reflected set `-x`.

### Input

- `x` -- interval

### Output

The `Interval` representing `-x`.

### Algorithm

If ``x = [a, b]``, then ``-x = [-b, -a]``.
"""
function reflect(x::Interval)
    return Interval(-max(x), -min(x))
end
