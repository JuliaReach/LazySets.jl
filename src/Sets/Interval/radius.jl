# in 1D, the radius is the same for any `p`
function radius(X::Interval{N}, ::Real=Inf) where {N}
    return diameter(X) / N(2)
end
