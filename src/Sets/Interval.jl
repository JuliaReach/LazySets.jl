using IntervalArithmetic: AbstractInterval, mince
import Base: +, -, *, ∈, ⊆, rand, min, max

export Interval,
       dim, σ, center,
       vertices_list,
       isflat,
       linear_map
       constraints_list

"""
    Interval{N, IN<:AbstractInterval{N}} <: AbstractHyperrectangle{N}

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
Interval{Float64, IntervalArithmetic.Interval{Float64}}([0, 1])
```
or a 2-vector:

```jldoctest interval_constructor
julia> x = Interval([0.0, 1.0])
Interval{Float64, IntervalArithmetic.Interval{Float64}}([0, 1])
```

Note that if the package `IntervalArithmetic` is loaded in the current scope,
you have to prepend `LazySets.` to the `Interval` type since there is a name
conflict otherwise.

```jldoctest interval_constructor
julia> using IntervalArithmetic
WARNING: using IntervalArithmetic.Interval in module Main conflicts with an existing identifier.

julia> x = Interval(IntervalArithmetic.Interval(0.0, 1.0))
Interval{Float64, IntervalArithmetic.Interval{Float64}}([0, 1])

julia> dim(x)
1

julia> center(x)
1-element Vector{Float64}:
 0.5
```

This type is such that the usual pairwise arithmetic operators `+`, `-`, `*` trigger
the corresponding interval arithmetic backend method, and return a new
`Interval` object. For the symbolic Minkowksi sum, use `MinkowskiSum` or `⊕`.

Interval of other numeric types can be created as well, eg. a rational interval:

```jldoctest interval_constructor
julia> Interval(0//1, 2//1)
Interval{Rational{Int64}, IntervalArithmetic.Interval{Rational{Int64}}}([0//1, 2//1])
```
"""
struct Interval{N, IN<:AbstractInterval{N}} <: AbstractHyperrectangle{N}
    dat::IN

    function Interval(dat::IN) where {N, IN<:AbstractInterval{N}}
        @assert isfinite(dat.lo) && isfinite(dat.hi) "intervals must be bounded"

        return new{N, IN}(dat)
    end
end

# constructor from two numbers with type promotion
function Interval(lo::N1, hi::N2) where {N1, N2}
    N = promote_type(N1, N2)
    Interval(IntervalArithmetic.Interval(N(lo), N(hi)))
end

# constructor from a vector
function Interval(x::AbstractVector)
    @assert length(x) == 2 "vector for Interval constructor has to be 2D"
    Interval(x[1], x[2])
end

# constructor from a single number
function Interval(x::Real)
    Interval(IntervalArithmetic.Interval(x))
end

isoperationtype(::Type{<:Interval}) = false
isconvextype(::Type{<:Interval}) = true

"""
    dim(x::Interval)

Return the ambient dimension of an interval.

### Input

- `x` -- interval

### Output

The integer 1.
"""
dim(x::Interval) = 1

"""
    σ(d::AbstractVector, x::Interval)

Return the support vector of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector, x::Interval)
    @assert length(d) == dim(x)
    N = promote_type(eltype(d), eltype(x))
    return d[1] > zero(N) ? high(x) : low(x)
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
    @assert length(d) == dim(x)
    N = promote_type(eltype(d), eltype(x))
    return d[1] * (d[1] > zero(N) ? max(x) : min(x))
end

"""
    center(x::Interval)

Return the interval's center.

### Input

- `x` -- interval

### Output

The center, or midpoint, of ``x``.
"""
function center(x::Interval)
    return [IntervalArithmetic.mid(x.dat)]
end

"""
    -(x::Interval, y::Interval)

Return the difference of the intervals.

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

Return the product of the intervals.

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

Return whether a vector is contained in the interval.

### Input

- `v` -- one-dimensional vector
- `x` -- interval

### Output

`true` iff `x` contains `v`'s first component.
"""
function ∈(v::AbstractVector, x::Interval)
    return v[1] ∈ x.dat
end

"""
    ∈(v::Number, x::Interval)

Return whether a number is contained in the interval.

### Input

- `v` -- scalar
- `x` -- interval

### Output

`true` iff `x` contains `v`.
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

The lower (`lo`) component of the interval.
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

The higher (`hi`) component of the interval.
"""
function max(x::Interval)
    return x.dat.hi
end

"""
    low(x::Interval{N}) where {N}

Return the lower coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the lower coordinate of the interval.
"""
function low(x::Interval{N}) where {N}
    return N[x.dat.lo]
end

"""
    high(x::Interval{N}) where {N}

Return the higher coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the higher coordinate of the interval.
"""
function high(x::Interval{N}) where {N}
    return N[x.dat.hi]
end

"""
    an_element(x::Interval)

Return some element of an interval.

### Input

- `x` -- interval

### Output

The left border (`min(x)`) of the interval.
"""
function an_element(x::Interval)
    return [min(x)]
end

"""
    rand(::Type{Interval}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
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
              seed::Union{Int, Nothing}=nothing)
    @assert dim == 1 "cannot create a random Interval of dimension $dim"
    rng = reseed(rng, seed)
    x = randn(rng, N)
    y = randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

"""
    vertices_list(x::Interval)

Return the list of vertices of this interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval represented as two one-dimensional vectors,
or just one if the interval is degenerate (the endpoints match within the tolerance).
"""
function vertices_list(x::Interval)
    a = min(x)
    b = max(x)
    return _isapprox(a, b) ? [[a]] : [[a], [b]]
end

"""
    constraints_list(x::Interval{N}) where {N}

Return the list of constraints of the given interval.

### Input

- `x` -- interval

### Output

The list of constraints of the interval represented as two one-dimensional
half-spaces.
"""
function constraints_list(x::Interval{N}) where {N}
    constraints = Vector{LinearConstraint{N, SingleEntryVector{N}}}(undef, 2)
    e₁ = SingleEntryVector(1, 1, one(N))
    constraints[1] = HalfSpace(e₁, max(x))
    constraints[2] = HalfSpace(-e₁, -min(x))
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

The in-place version is not available because interval types
are immutable.
"""
function translate(x::Interval, v::AbstractVector)
    @assert length(v) == dim(x) "cannot translate a $(dim(x))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Interval(x.dat + v[1])
end


# --- AbstractHyperrectangle interface functions ---


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
    @assert i == 1 "an interval is one-dimensional"
    return (max(x) - min(x)) / N(2)
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
    return [radius_hyperrectangle(x, 1)]
end

"""
    isflat(I::Interval)

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
function isflat(I::Interval)
    return isapproxzero(IntervalArithmetic.diam(I.dat))
end

"""
    plot_recipe(I::Interval{N}, [ε]=zero(N)) where {N}

Convert an interval to a pair `(x, y)` of points for plotting.

### Input

- `I` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.

### Notes

We consider the interval as a line segment with y coordinate equal to zero.
"""
function plot_recipe(I::Interval{N}, ε=zero(N)) where {N}
    return [min(I), max(I)], zeros(N, 2)
end

"""
    linear_map(M::AbstractMatrix, x::Interval)

Concrete linear map of an interval.

### Input

- `M` -- matrix
- `x` -- interval

### Output

Either an interval or a zonotope, depending on the leading dimension (i.e. the
number of rows) of `M`:

- If `size(M, 1) == 1`, the output is an interval obtained by scaling `x` by the
  matrix `M`.
- If `size(M, 1) > 1`, the output is a zonotope whose center is `M * center(x)`
  and it has the single generator, `M * g`, where `g = (high(x)-low(x))/2`.
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
    α = M[1, 1]
    return Interval(α * x.dat)
end

function _linear_map_zonotope(M::AbstractMatrix, x::Interval)
    nout = size(M, 1)
    cx = IntervalArithmetic.mid(x.dat)
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

The interval obtained by applying the numerical scale to the given interval.
"""
function scale(α::Real, x::Interval)
    return Interval(α * x.dat)
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
    @boundscheck i == 1 || throw(ArgumentError("an interval has dimension one, but the index is $i"))
    return IntervalArithmetic.mid(x.dat)
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

Compute the diameter of an interval, defined as ``\\Vert b - a\\Vert`` in the ``p`-norm, where ``a`` (resp. ``b``) are the minimum (resp. maximum) of the given interval.

### Input

- `x` -- interval
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.

### Notes

In one dimension all norms are the same.
"""
function diameter(x::Interval, p::Real=Inf)
    return IntervalArithmetic.diam(x.dat)
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
    return split(x, k[1])
end

function split(x::Interval, k::Int)
    @assert k > 0 "can only split into a positive number of intervals"
    return [Interval(x2) for x2 in mince(x.dat, k)]
end

function vertices_list(x::IntervalArithmetic.Interval{N}) where {N}
    a = IntervalArithmetic.inf(x)
    b = IntervalArithmetic.sup(x)
    ST = IntervalArithmetic.SVector{1, N}
    if _isapprox(a, b)
        vlist = [ST(a)]
    else
        vlist = [ST(a), ST(b)]
    end
    return vlist
end

function vertices_list(H::IntervalArithmetic.IntervalBox)
    return vertices_list(convert(Hyperrectangle, H))
end

function chebyshev_center(x::Interval{N}; compute_radius::Bool=false) where {N}
    if compute_radius
        return center(x), zero(N)
    end
    return center(x)
end

"""
    fast_interval_pow(a::IA.Interval, n::Int)

Computes the `n`th power of the given interval `a` without using correct rounding.

### Input

- `a` -- interval (from `IntervalArithmetic.jl`)
- `n` -- integer

### Output

A non-rigorous approximation of `a^n`.

### Notes

For a rigorous approximation with correct rounding,
use `a^n` from `IntervalArithmetic.jl`.

Review after IntervalArithmetic.jl#388
"""
const IA = IntervalArithmetic
function fast_interval_pow(a::IA.Interval, n::Int)
    if iszero(n)
        return one(a)
    elseif isodd(n)
        return IA.Interval(a.lo ^ n, a.hi ^ n)
    else
        if 0 ∈ a
            return IA.Interval(zero(a.lo), max(abs(a.lo), abs(a.hi)) ^ n)
        else
            lon = a.lo ^ n
            hin = a.hi ^ n
            return IA.Interval(min(lon, hin), max(lon, hin))
        end
    end
end
