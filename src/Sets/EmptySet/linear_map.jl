@validate function linear_map(M::AbstractMatrix, ∅::EmptySet)
    N = promote_type(eltype(M), eltype(∅))
    return EmptySet{N}(size(M, 1))
end
