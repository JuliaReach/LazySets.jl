function linear_map_inverse(Minv::AbstractMatrix{N}, ∅::EmptySet{N}) where {N}
    @assert size(Minv, 1) == dim(∅) "a linear map of size $(size(Minv)) " *
                                    "cannot be applied to a set of dimension $(dim(∅))"
    n = size(Minv, 2)
    return EmptySet{N}(n)
end
