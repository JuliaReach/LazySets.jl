"""
    convert(::Type{Interval}, X::LazySet)

Convert a one-dimensional convex set to an interval.

### Input

- `Interval` -- interval target type
- `X`        -- one-dimensional convex set

### Output

An interval.

### Examples

```jldoctest
julia> convert(Interval, Hyperrectangle([0.5], [0.5]))
Interval{Float64}([0, 1])
```
"""
function convert(::Type{Interval}, X::LazySet)
    return Interval(convert(IA.Interval, X))
end

"""
    convert(::Type{IntervalArithmetic.Interval}, X::LazySet)

Convert a convex set to an `Interval` from `IntervalArithmetic`.

### Input

- `Interval` -- target type, from `IntervalArithmetic`
- `X`        -- convex set

### Output

An `IntervalArithmetic.Interval`.
"""
function convert(::Type{IA.Interval}, X::LazySet)
    @assert dim(X) == 1 "cannot convert a $(dim(X))-dimensional set to an `Interval`"
    @assert isconvextype(typeof(X)) "cannot convert a non-convex set to an `Interval`"

    l, h = extrema(X, 1)
    return IA.interval(l, h)
end

"""
    convert(::Type{Interval}, x::IntervalArithmetic.Interval)

Convert an `Interval` from `IntervalArithmetic` to an `Interval` in `LazySets`.

### Input

- `Interval` -- target type
- `x`        -- interval (`IntervalArithmetic.Interval`)

### Output

A `LazySets.Interval`.
"""
function convert(::Type{Interval}, x::IA.Interval)
    return Interval(x)
end
