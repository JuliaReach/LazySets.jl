function center(X::Interval)
    return [_center(X)]
end

@validate function center(X::Interval, i::Int)
    return _center(X)
end

# returns a number, not a vector
function _center(X::Interval)
    return IA.mid(X.dat)
end
