function convex_hull(U1::Universe, U2::Universe)
    return _convex_hull_universe(U1, U2)
end

function _convex_hull_universe(U::Universe, X::LazySet)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"
    return U
end
