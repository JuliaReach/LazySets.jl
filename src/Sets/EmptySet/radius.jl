function radius(∅::EmptySet, ::Real=Inf)
    N = eltype(∅)
    return zero(N)
end
