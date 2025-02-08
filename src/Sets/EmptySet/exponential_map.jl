function exponential_map(M::AbstractMatrix, ∅::EmptySet)
    n = dim(∅)
    @assert size(M) == (n, n) "cannot apply an exponential map of dimension " *
                              "$(size(M)) to an $n-dimensional set"

    N = promote_type(eltype(M), eltype(∅))
    return EmptySet{N}(size(M, 1))
end
