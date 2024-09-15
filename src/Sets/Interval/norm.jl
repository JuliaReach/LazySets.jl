# in 1D, the norm is the same for any `p`
function norm(X::Interval, ::Real=Inf)
    return max(abs(min(X)), abs(max(X)))
end
