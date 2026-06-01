@validate function isequivalent(U1::Universe, U2::Universe, witness::Bool=false)
    return _witness_result_empty(witness, true, U1, U2)
end
