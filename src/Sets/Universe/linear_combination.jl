@validate function linear_combination(U1::Universe, U2::Universe)
    return _linear_combination_universe(U1, U2)
end

function _linear_combination_universe(U::Universe, X::LazySet)
    if isempty(X)
        return X
    end
    return U
end
