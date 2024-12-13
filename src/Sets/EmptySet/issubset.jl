function ⊆(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _issubset_emptyset(∅₁, ∅₂, witness)
end

function _issubset_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return _witness_result_empty(witness, true, ∅, X)
end

function _issubset_emptyset2(X::LazySet, ∅::EmptySet, witness::Bool=false)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"

    result = isempty(X)
    if !result && witness
        return _witness_result_empty(witness, false, X, ∅, an_element(X))
    else
        return _witness_result_empty(witness, result, X, ∅)
    end
end
