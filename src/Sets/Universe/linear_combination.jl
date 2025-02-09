function linear_combination(U1::Universe, U2::Universe)
    return _linear_combination_universe(U1, U2)
end

function _linear_combination_universe(U::Universe, X::LazySet)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"

    if isempty(X)
        return X
    end
    return U
end
