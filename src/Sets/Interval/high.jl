function high(X::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return max(X)
end

function high(X::Interval)
    return [max(X)]
end
