@validate function radius(∅::EmptySet, p::Real=Inf)
    N = eltype(∅)
    return zero(N)
end
