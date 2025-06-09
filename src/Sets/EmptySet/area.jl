function area(∅::EmptySet)
    @assert dim(∅) ∈ (2, 3) "this function only applies to two-dimensional " *
                              "or three-dimensional sets, but the given set " *
                              "is $(dim(∅))-dimensional"

    N = eltype(∅)
    return zero(N)
end
