"""
    linear_map(M::AbstractMatrix{N}, ∅::EmptySet{N}) where {N}

Return the linear map of an empty set.

### Input

- `M` -- matrix
- `∅` -- empty set

### Output

An empty set.
"""
function linear_map(M::AbstractMatrix, ∅::EmptySet)
    N = eltype(∅)
    @assert size(M, 2) == dim(∅) "cannot apply a $(size(M))-dimensional " *
                                 "matrix to a $(dim(∅))-dimensional set"

    return EmptySet{N}(size(M, 1))
end
