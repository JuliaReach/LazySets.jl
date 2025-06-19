function center(X::Interval)
    return [_center(X)]
end

@inline function center(X::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return _center(X)
end

# returns a number, not a vector
function _center(X::Interval)
    return IA.mid(X.dat)
end
