@validate function difference(∅₁::EmptySet, ∅₂::EmptySet)
    return _difference_emptyset(∅₁, ∅₂)
end

function _difference_emptyset(∅::EmptySet, ::LazySet)
    return ∅
end

function _difference_emptyset2(X::LazySet, ::EmptySet)
    return X
end
