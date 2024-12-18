function linear_map(M::AbstractMatrix, Z::ZeroSet)
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"

    N = promote_type(eltype(M), eltype(Z))
    return ZeroSet{N}(size(M, 1))
end
