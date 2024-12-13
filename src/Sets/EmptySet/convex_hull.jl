function convex_hull(∅::EmptySet)
    return ∅
end

function convex_hull(∅₁::EmptySet, ∅₂::EmptySet)
    return _convex_hull_emptyset(∅₁, ∅₂)
end

function _convex_hull_emptyset(∅::EmptySet, X::EmptySet)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return ∅
end
