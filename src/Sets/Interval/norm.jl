function norm(X::Interval, p::Real=Inf)
    return max(abs(min(X)), abs(max(X)))
end
