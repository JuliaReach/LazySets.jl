@validate function convex_hull(I1::Interval, I2::Interval)
    return Interval(min(_min(I1), _min(I2)), max(_max(I1), _max(I2)))
end
