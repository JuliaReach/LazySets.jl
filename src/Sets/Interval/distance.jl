function distance(X::Interval, Y::Interval; p::Real=2)
    d = max(min(X) - max(Y), min(Y) - max(X))
    if d < 0
        return zero(d)
    end
    return d
end
