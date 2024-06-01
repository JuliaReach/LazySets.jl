function ≈(I1::Interval, I2::Interval)
    return _isapprox(min(I1), min(I2)) && _isapprox(max(I1), max(I2))
end
