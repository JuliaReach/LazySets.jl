@validate function convex_hull(I1::Interval, I2::Interval)
    return Interval(min(min(I1), min(I2)), max(max(I1), max(I2)))
end
