function linear_map(M::AbstractMatrix, S::Singleton)
    @assert dim(S) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(S))"

    return Singleton(M * S.element)
end
