@validate function minkowski_sum(U1::Universe, U2::Universe)
    return _minkowski_sum_universe(U1, U2)
end

function _minkowski_sum_universe(U::Universe, X::LazySet)
    if isempty(X)
        return X
    end
    return U
end
