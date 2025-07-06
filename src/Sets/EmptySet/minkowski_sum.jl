@validate function minkowski_sum(∅₁::EmptySet, ∅₂::EmptySet)
    return _minkowski_sum_emptyset(∅₁, ∅₂)
end

function _minkowski_sum_emptyset(∅::EmptySet, X::LazySet)
    return ∅
end
