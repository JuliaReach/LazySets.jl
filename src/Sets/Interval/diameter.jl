# in one dimension, all p-norms are identical
function diameter(X::Interval, ::Real=Inf)
    return max(X) - min(X)
end
