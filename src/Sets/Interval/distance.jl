@validate function distance(X::Interval, Y::Interval; p::Real=2)
    d = max(_min(X) - _max(Y), _min(Y) - _max(X))
    if d < 0
        return zero(d)
    end
    return d
end
