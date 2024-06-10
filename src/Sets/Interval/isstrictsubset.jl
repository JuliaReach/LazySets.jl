function ⊂(X::Interval, Y::Interval, witness::Bool=false)
    if min(X) < min(Y) || max(X) > max(Y)
        return _witness_result_empty(witness, false, X, Y)
    end
    if min(X) > min(Y)
        return witness ? (true, low(Y)) : true
    elseif max(Y) > max(X)
        return witness ? (true, high(Y)) : true
    end
    return _witness_result_empty(witness, false, X, Y)
end
