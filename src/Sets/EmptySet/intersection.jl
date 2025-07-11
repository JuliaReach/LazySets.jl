@validate function intersection(∅₁::EmptySet, ∅₂::EmptySet)
    return _intersection_emptyset(∅₁::EmptySet, ∅₂::LazySet)
end

function _intersection_emptyset(∅::EmptySet, X::LazySet)
    return ∅
end
