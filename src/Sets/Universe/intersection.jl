function intersection(U1::Universe, U2::Universe)
    return _intersection_universe(U1, U2)
end

function _intersection_universe(U::Universe, X::LazySet)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"
    return X
end
