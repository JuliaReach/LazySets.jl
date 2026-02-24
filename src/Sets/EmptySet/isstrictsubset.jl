@validate function ⊂(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    return _isstrictsubset_emptyset(∅₁, ∅₂, witness)
end

function _isstrictsubset_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    if isempty(X)
        return _witness_result_empty(witness, false, ∅, X)
    else
        witness ? (true, an_element(X)) : true
    end
end

function _isstrictsubset_emptyset2(X::LazySet, ∅::EmptySet, witness::Bool=false)
    if witness
        empty, w = isempty(X, true)
        if empty
            N = promote_type(eltype(X), eltype(∅))
            return (false, N[])
        else
            return _witness_result(witness, false, w)
        end
    end
    return false
end
