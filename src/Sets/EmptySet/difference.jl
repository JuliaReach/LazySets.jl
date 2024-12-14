function difference(∅₁::EmptySet, ∅₂::EmptySet)
    return _difference_emptyset(∅₁, ∅₂)
end

function _difference_emptyset(∅::EmptySet, X::LazySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return ∅
end

function _difference_emptyset2(X::LazySet, ∅::EmptySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return X
end
