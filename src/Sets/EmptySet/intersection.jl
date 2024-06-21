function intersection(∅₁::EmptySet, ∅₂::EmptySet)
    @assert dim(∅₁) == dim(∅₂) "cannot take the intersection between two " *
                               "empty sets of dimensions $(dim(∅₁)) and $(dim(∅₂))"
    return ∅₁
end
