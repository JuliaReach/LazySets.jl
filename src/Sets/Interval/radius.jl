# in 1D, the radius is the same for any `p`
@validate function radius(X::Interval{N}, p::Real=Inf) where {N}
    return diameter(X) / N(2)
end
