# in 1D, the norm is the same for any `p`
@validate function norm(X::Interval, p::Real=Inf)
    return max(abs(min(X)), abs(max(X)))
end
