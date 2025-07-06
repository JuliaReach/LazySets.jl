function scale(α::Real, X::Interval)
    return Interval(α * X.dat)
end
