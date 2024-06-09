# the implementations here are equivalent to the default implementation for
# `LazySet`, but more efficient than the more specific implementations for
# subtypes

function extrema(X::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return (min(X), max(X))
end

function extrema(X::Interval)
    return ([min(X)], [max(X)])
end
