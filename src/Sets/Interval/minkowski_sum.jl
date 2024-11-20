function minkowski_sum(x::Interval, y::Interval)
    return Interval(x.dat + y.dat)
end
