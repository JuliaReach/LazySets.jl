@validate function minkowski_difference(∅₁::EmptySet, ∅₂::EmptySet)
    return _minkowski_difference_emptyset(∅₁, ∅₂)
end

function _minkowski_difference_emptyset(∅::EmptySet, X::LazySet)
    return ∅
end

function _minkowski_difference_emptyset2(X::LazySet, ∅::EmptySet)
    return X
end
