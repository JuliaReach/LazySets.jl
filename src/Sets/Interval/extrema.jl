# the implementations here are equivalent to the default implementation for
# `LazySet`, but more efficient than the more specific implementations for
# subtypes

@validate function extrema(X::Interval, i::Int)
    return (_min(X), _max(X))
end

function extrema(X::Interval)
    return ([_min(X)], [_max(X)])
end
