isdisjoint(∅1::EmptySet, ∅2::EmptySet, witness::Bool=false) = _isdisjoint_emptyset(∅1, ∅2, witness)

function _isdisjoint_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return _witness_result_empty(witness, true, ∅, X)
end
