function ⊆(U1::Universe, U2::Universe, witness::Bool=false)
    return _issubset_universe2(U1, U2, witness)
end

function _issubset_universe(U::Universe, X::LazySet, witness::Bool=false)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"

    return isuniversal(X, witness)
end

function _issubset_universe2(X::LazySet, U::Universe, witness::Bool=false)
    @assert dim(X) == dim(U) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(U)), respectively"

    return _witness_result_empty(witness, true, X, U)
end
