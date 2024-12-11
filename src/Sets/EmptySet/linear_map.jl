function linear_map(M::AbstractMatrix, ∅::EmptySet)
    @assert size(M, 2) == dim(∅) "cannot apply a $(size(M))-dimensional " *
                                 "matrix to a $(dim(∅))-dimensional set"

    N = eltype(∅)
    return EmptySet{N}(size(M, 1))
end
