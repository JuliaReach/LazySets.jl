function convex_hull(∅::EmptySet)
    return ∅
end

@validate function convex_hull(∅₁::EmptySet, ∅₂::EmptySet)
    return _convex_hull_emptyset(∅₁, ∅₂)
end

function _convex_hull_emptyset(∅::EmptySet, X::LazySet)
    return convex_hull(X)
end
