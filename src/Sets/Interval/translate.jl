"""
# Extended help

    translate(X::Interval, v::AbstractVector)

### Notes

An in-place version is not available because the `IntervalArithmetic.Interval`
type is immutable.
"""
@validate function translate(X::Interval, v::AbstractVector)
    return Interval(X.dat + @inbounds v[1])
end
