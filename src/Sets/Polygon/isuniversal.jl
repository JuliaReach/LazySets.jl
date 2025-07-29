function isuniversal(P::Polygon, witness::Bool=false)
    return witness ? (false, _non_element(P)) : false
end

function _non_element(P::Polygon)
    return high(P) + ones(eltype(P), 2)
end
