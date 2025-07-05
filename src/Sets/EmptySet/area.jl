@validate function area(∅::EmptySet)
    N = eltype(∅)
    return zero(N)
end
