# in one dimension, all p-norms are identical
@validate function diameter(X::Interval, p::Real=Inf)
    return _max(X) - _min(X)
end
