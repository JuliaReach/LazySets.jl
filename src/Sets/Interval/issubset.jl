function ⊆(x::Interval, y::Interval, witness::Bool=false)
    if min(y) > min(x)
        return witness ? (false, low(x)) : false
    elseif max(x) > max(y)
        return witness ? (false, high(x)) : false
    end
    return _witness_result_empty(witness, true, x, y)
end
