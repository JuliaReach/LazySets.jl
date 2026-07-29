function isapprox(I1::Interval, I2::Interval)
    return _isapprox(_min(I1), _min(I2)) && _isapprox(_max(I1), _max(I2))
end
