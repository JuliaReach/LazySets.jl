function cartesian_product(∅₁::EmptySet, ∅₂::EmptySet)
    return _cartesian_product_emptyset(∅₁, ∅₂)
end

function _cartesian_product_emptyset(∅::EmptySet, X::LazySet)
    N = promote_type(eltype(∅), eltype(X))
    return EmptySet{N}(dim(∅) + dim(X))
end
