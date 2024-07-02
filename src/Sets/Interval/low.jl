function low(X::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return min(X)
end

function low(X::Interval)
    return [min(X)]
end
