function ⊆(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _issubset_emptyset(∅₁, ∅₂, witness)
end

function _issubset_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return _witness_result_empty(witness, true, ∅, X)
end
