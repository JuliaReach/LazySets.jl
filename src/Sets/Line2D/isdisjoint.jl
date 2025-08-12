@validate function isdisjoint(L1::Line2D, L2::Line2D, witness::Bool=false)
    disjoint = _isdisjoint(L1, L2)
    if disjoint
        return _witness_result_empty(witness, true, L1, L2)
    end
    return witness ? (false, an_element(intersection(L1, L2))) : false
end

# the lines do not intersect <=> det is zero and they are not identical
function _isdisjoint(L1::Line2D, L2::Line2D)
    det = right_turn(L1.a, L2.a)
    disjoint = isapproxzero(det) && !_isapprox(L1.b, L2.b)
    return disjoint
end
