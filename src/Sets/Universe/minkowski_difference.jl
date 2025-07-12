@validate function minkowski_difference(U1::Universe, U2::Universe)
    return _minkowski_difference_universe(U1, U2)
end

function _minkowski_difference_universe(U::Universe, X::LazySet)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"
    return U
end

function _minkowski_difference_universe2(X::LazySet, U::Universe)
    @assert dim(X) == dim(U) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(U)), respectively"

    if isuniversal(X)
        return U
    end
    N = promote_type(eltype(X), eltype(U))
    return EmptySet{N}(dim(U))
end
