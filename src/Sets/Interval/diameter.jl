# in one dimension, all p-norms are identical
function diameter(x::Interval, ::Real=Inf)
    return max(x) - min(x)
end
