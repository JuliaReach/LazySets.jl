@validate function isdisjoint(U1::Universe, U2::Universe, witness::Bool=false)
    return _isdisjoint_universe(U1, U2, witness)
end

function _isdisjoint_universe(U::Universe, X::LazySet, witness)
    if isempty(X)
        return _witness_result_empty(witness, true, U, X)
    else
        return witness ? (false, an_element(X)) : false
    end
end
