@validate function linear_combination(∅₁::EmptySet, ∅₂::EmptySet)
    return _linear_combination_emptyset(∅₁, ∅₂)
end

# the empty set is absorbing for the linear combination
function _linear_combination_emptyset(∅::EmptySet, ::LazySet)
    return ∅
end
