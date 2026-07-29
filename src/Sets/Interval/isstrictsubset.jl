@validate function ⊂(X::Interval, Y::Interval, witness::Bool=false)
    if _min(X) < _min(Y)
        return witness ? (false, low(X)) : false
    elseif _max(X) > _max(Y)
        return witness ? (false, high(X)) : false
    end
    if _min(X) > _min(Y)
        return witness ? (true, low(Y)) : true
    elseif _max(Y) > _max(X)
        return witness ? (true, high(Y)) : true
    end
    return _witness_result_empty(witness, false, X, Y)
end
