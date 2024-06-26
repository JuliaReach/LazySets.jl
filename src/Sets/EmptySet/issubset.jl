function ⊆(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _witness_result_empty(witness, true, ∅₁, ∅₂)
end
