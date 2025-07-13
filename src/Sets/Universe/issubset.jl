@validate function ⊆(U1::Universe, U2::Universe, witness::Bool=false)
    return _issubset_universe2(U1, U2, witness)
end

function _issubset_universe(U::Universe, X::LazySet, witness::Bool=false)
    return isuniversal(X, witness)
end

function _issubset_universe2(X::LazySet, U::Universe, witness::Bool=false)
    return _witness_result_empty(witness, true, X, U)
end
