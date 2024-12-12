function cartesian_product(∅₁::EmptySet, ∅₂::EmptySet)
    return _cartesian_product_emptyset(∅₁, ∅₂)
end

function _cartesian_product_emptyset(∅::EmptySet, X::LazySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    throw(ArgumentError("cannot take the Cartesian product with an empty set"))
end
