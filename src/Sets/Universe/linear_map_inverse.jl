function linear_map_inverse(Minv::AbstractMatrix{N}, U::Universe{N}) where {N}
    @assert size(Minv, 1) == dim(U) "a linear map of size $(size(Minv)) " *
                                    "cannot be applied to a universe of dimension $(dim(U))"
    n = size(Minv, 2)
    return Universe{N}(n)
end
