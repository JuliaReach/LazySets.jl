@validate function norm(∅::EmptySet, p::Real=Inf)
    N = eltype(∅)
    return zero(N)
end
