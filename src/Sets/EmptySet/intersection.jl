function intersection(∅₁::EmptySet, ∅₂::EmptySet)
    return _intersection_emptyset(∅₁::EmptySet, ∅₂::LazySet)
end

function _intersection_emptyset(∅::EmptySet, X::LazySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return ∅
end
