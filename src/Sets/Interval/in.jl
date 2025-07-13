@validate function ∈(v::AbstractVector, X::Interval)
    return @inbounds v[1] ∈ X.dat
end
