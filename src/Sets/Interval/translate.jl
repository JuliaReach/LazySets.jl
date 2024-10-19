"""
# Extended help

    translate(x::Interval, v::AbstractVector)

### Notes

An in-place version is not available because the `IntervalArithmetic.Interval`
type is immutable.
"""
function translate(x::Interval, v::AbstractVector)
    @assert length(v) == dim(x) "cannot translate a $(dim(x))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Interval(x.dat + @inbounds v[1])
end
