@validate function minkowski_sum(X::Interval, Y::Interval)
    return Interval(X.dat + Y.dat)
end
