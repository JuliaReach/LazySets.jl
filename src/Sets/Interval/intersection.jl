"""
    intersection(x::Interval, y::Interval)

Compute the intersection of two intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

If the intervals do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection.
"""
function intersection(x::Interval, y::Interval)
    l = max(min(x), min(y))
    h = min(max(x), max(y))
    if l > h
        require(@__MODULE__, :LazySets; fun_name="intersection")

        N = promote_type(eltype(x), eltype(y))
        return EmptySet{N}(1)
    else
        return Interval(l, h)
    end
end
