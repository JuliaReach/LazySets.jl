@validate function intersection(U1::Universe, U2::Universe)
    return _intersection_universe(U1, U2)
end

function _intersection_universe(U::Universe, X::LazySet)
    return X
end
