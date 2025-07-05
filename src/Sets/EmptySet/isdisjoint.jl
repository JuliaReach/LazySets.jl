@validate function isdisjoint(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _isdisjoint_emptyset(∅₁, ∅₂, witness)
end

function _isdisjoint_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _witness_result_empty(witness, true, ∅, X)
end
