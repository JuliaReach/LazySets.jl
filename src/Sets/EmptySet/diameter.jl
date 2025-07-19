@validate function diameter(∅::EmptySet, p::Real=Inf)
    N = eltype(∅)
    return zero(N)
end
