@validate function high(X::Interval, i::Int)
    return _max(X)
end

function high(X::Interval)
    return [_max(X)]
end
