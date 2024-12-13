function cartesian_product(∅₁::EmptySet, ∅₂::EmptySet)
    return _cartesian_product_emptyset(∅₁, ∅₂)
end

function _cartesian_product_emptyset(∅::EmptySet, X::LazySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    N = promote_type(eltype(∅), eltype(X))
    return EmptySet{N}(dim(∅) + dim(X))
end
