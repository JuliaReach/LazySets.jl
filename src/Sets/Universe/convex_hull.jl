@validate function convex_hull(U1::Universe, U2::Universe)
    return _convex_hull_universe(U1, U2)
end

function _convex_hull_universe(U::Universe, X::LazySet)
    return U
end
