function minkowski_sum(∅₁::EmptySet, ∅₂::EmptySet)
    return _minkowski_sum_emptyset(∅₁, ∅₂)
end

function _minkowski_sum_emptyset(∅::EmptySet, X::LazySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return ∅
end
