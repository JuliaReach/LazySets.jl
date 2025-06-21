function isempty(∅::EmptySet, witness::Bool=false)
    return _witness_result_empty(witness, true, eltype(∅))
end
