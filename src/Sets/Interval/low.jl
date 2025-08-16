@validate function low(X::Interval, i::Int)
    return min(X)
end

function low(X::Interval)
    return [min(X)]
end
