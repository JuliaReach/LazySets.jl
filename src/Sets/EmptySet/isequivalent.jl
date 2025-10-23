@validate function isequivalent(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _issubset_emptyset2(∅₁, ∅₂, witness)
end
