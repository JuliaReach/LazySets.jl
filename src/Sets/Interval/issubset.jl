@validate function ⊆(X::Interval, Y::Interval, witness::Bool=false)
    if min(Y) > min(X)
        return witness ? (false, low(X)) : false
    elseif max(X) > max(Y)
        return witness ? (false, high(X)) : false
    end
    return _witness_result_empty(witness, true, X, Y)
end
