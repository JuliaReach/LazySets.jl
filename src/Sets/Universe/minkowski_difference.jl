@validate function minkowski_difference(U1::Universe, U2::Universe)
    return _minkowski_difference_universe(U1, U2)
end

function _minkowski_difference_universe(U::Universe, X::LazySet)
    return U
end

function _minkowski_difference_universe2(X::LazySet, U::Universe)
    if isuniversal(X)
        return U
    end
    N = promote_type(eltype(X), eltype(U))
    return EmptySet{N}(dim(U))
end
