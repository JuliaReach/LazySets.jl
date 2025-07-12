@validate function difference(U1::Universe, U2::Universe)
    return _difference_universe2(U1, U2)
end

function _difference_universe(U::Universe, X::LazySet)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"
    return complement(X)
end

function _difference_universe2(X::LazySet, U::Universe)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"

    N = promote_type(eltype(X), eltype(U))
    return EmptySet{N}(dim(U))
end
