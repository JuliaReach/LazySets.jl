function intersection(x::Interval, y::Interval)
    l = max(min(x), min(y))
    h = min(max(x), max(y))
    if l > h
        N = promote_type(eltype(x), eltype(y))
        return EmptySet{N}(1)
    else
        return Interval(l, h)
    end
end
