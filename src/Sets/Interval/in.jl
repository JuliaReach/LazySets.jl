@validate function in(v::AbstractVector, X::Interval)
    return @inbounds v[1] ∈ X.dat
end
