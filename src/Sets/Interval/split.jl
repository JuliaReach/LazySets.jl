"""
    split(X::Interval, k)

Partition an interval into `k` uniform sub-intervals.

### Input

- `X` -- interval
- `k` -- number of sub-intervals, possibly wrapped in a vector of length 1

### Output

A list of `k` `Interval`s.
"""
function split(X::Interval, k::AbstractVector{Int})
    @assert length(k) == 1 "an interval can only be split along one dimension"
    return split(X, @inbounds k[1])
end

function split(X::Interval, k::Int)
    @assert k > 0 "can only split into a positive number of intervals"
    return [Interval(x2) for x2 in IA.mince(X.dat, k)]
end
