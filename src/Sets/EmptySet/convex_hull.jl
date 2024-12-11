function convex_hull(∅::EmptySet)
    return ∅
end

function convex_hull(∅₁::EmptySet, ∅₂::EmptySet)
    @assert dim(∅₁) == dim(∅₂) "cannot take the convex hull between two " *
                               "empty sets of dimensions $(dim(∅₁)) and $(dim(∅₂))"
    return ∅₁
end
