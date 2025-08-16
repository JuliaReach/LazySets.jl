@validate function high(X::Interval, i::Int)
    return max(X)
end

function high(X::Interval)
    return [max(X)]
end
