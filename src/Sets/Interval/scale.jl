function scale(α::Real, x::Interval)
    return Interval(α * x.dat)
end
