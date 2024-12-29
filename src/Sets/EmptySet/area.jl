function area(∅::EmptySet)
    @assert dim(∅) == 2 "this function only applies to two-dimensional sets, " *
                        "but the given set is $(dim(∅))-dimensional"

    N = eltype(∅)
    return zero(N)
end
