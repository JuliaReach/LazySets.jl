@validate function in(v::AbstractVector, X::Interval)
    e = @inbounds v[1]
    return IA.in_interval(e, X.dat)
end
