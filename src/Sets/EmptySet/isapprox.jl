function ≈(∅₁::EmptySet, ∅₂::EmptySet)
    @assert dim(∅₁) == dim(∅₂) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅₁)) and $(dim(∅₂)), respectively"
    return true
end
