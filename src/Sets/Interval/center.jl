function center(x::Interval)
    return [_center(x)]
end

@inline function center(x::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return _center(x)
end

# returns a number, not a vector
function _center(x::Interval)
    return IA.mid(x.dat)
end
