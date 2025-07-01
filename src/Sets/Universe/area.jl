function area(U::Universe)
    n = dim(U)
    @assert n ∈ (2, 3) "this function only applies to two-dimensional or " *
                        "three-dimensional sets, but the given set is " *
                        "$n-dimensional"

    N = eltype(U)
    return N(Inf)
end
