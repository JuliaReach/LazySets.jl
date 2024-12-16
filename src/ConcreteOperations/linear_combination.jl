function linear_combination(X::ConvexSet, Y::ConvexSet)
    return convex_hull(X, Y)
end

@commutative function linear_combination(∅::EmptySet, X::LazySet)
    return _linear_combination_emptyset(∅, X)
end

# disambiguation
@commutative function linear_combination(∅::EmptySet, X::ConvexSet)
    return _linear_combination_emptyset(∅, X)
end
