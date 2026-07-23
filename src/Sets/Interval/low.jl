@validate function low(X::Interval, i::Int)
    return _min(X)
end

function low(X::Interval)
    return [_min(X)]
end
