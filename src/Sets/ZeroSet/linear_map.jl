@validate function linear_map(M::AbstractMatrix, Z::ZeroSet)
    N = promote_type(eltype(M), eltype(Z))
    return ZeroSet{N}(size(M, 1))
end
