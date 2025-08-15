# in one dimension, all p-norms are identical
@validate function diameter(X::Interval, p::Real=Inf)
    return max(X) - min(X)
end
