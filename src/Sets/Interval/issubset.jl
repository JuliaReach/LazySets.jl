@validate function issubset(X::Interval, Y::Interval, witness::Bool=false)
    if _min(Y) > _min(X)
        return witness ? (false, low(X)) : false
    elseif _max(X) > _max(Y)
        return witness ? (false, high(X)) : false
    end
    return _witness_result_empty(witness, true, X, Y)
end
