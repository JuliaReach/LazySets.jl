@validate function issubset(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _issubset_emptyset(∅₁, ∅₂, witness)
end

function _issubset_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _witness_result_empty(witness, true, ∅, X)
end

function _issubset_emptyset2(X::LazySet, ∅::EmptySet, witness::Bool=false)
    if isempty(X)
        return _witness_result_empty(witness, true, X, ∅)
    else
        return witness ? (false, an_element(X)) : false
    end
end
