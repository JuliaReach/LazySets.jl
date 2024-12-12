function linear_map(M::AbstractMatrix, ∅::EmptySet)
    @assert size(M, 2) == dim(∅) "cannot apply a $(size(M))-dimensional " *
                                 "matrix to a $(dim(∅))-dimensional set"

    N = promote_type(eltype(M), eltype(∅))
    return EmptySet{N}(size(M, 1))
end
