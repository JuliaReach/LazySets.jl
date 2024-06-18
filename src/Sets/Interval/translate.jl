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
