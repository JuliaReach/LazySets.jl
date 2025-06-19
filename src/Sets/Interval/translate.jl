"""
# Extended help

    translate(X::Interval, v::AbstractVector)

### Notes

An in-place version is not available because the `IntervalArithmetic.Interval`
type is immutable.
"""
function translate(X::Interval, v::AbstractVector)
    @assert length(v) == dim(X) "cannot translate a $(dim(X))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Interval(X.dat + @inbounds v[1])
end
