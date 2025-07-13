@validate function area(U::Universe)
    N = eltype(U)
    return N(Inf)
end
