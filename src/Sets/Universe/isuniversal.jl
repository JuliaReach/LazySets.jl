function isuniversal(U::Universe, witness::Bool=false)
    return _witness_result_empty(witness, true, eltype(U))
end
