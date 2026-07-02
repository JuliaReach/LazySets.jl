@validate function minkowski_difference(U1::Universe, U2::Universe)
    return _minkowski_difference_universe(U1, U2)
end

function _minkowski_difference_universe(U::Universe, X::LazySet)
    return U
end

# see ext/LazySets/LazySetsUniverseExt.jl
_minkowski_difference_universe2(X, U) = error
