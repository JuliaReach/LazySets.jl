function linear_combination(X::LazySet, Y::LazySet)
    if isconvextype(typeof(X)) && isconvextype(typeof(Y))
        return _linear_combination_convex(X, Y)
    end
    throw(ArgumentError("the linear combination of non-convex sets is not implemented"))
end

function linear_combination(X::ConvexSet, Y::ConvexSet)
    return _linear_combination_convex(X, Y)
end

function _linear_combination_convex(X, Y)
    @assert dim(X) == dim(Y) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(Y)), respectively"

    if isempty(X)
        return X
    elseif isempty(Y)
        return Y
    end
    return convex_hull(X, Y)
end

@commutative function linear_combination(∅::EmptySet, X::LazySet)
    return _linear_combination_emptyset(∅, X)
end

@commutative function linear_combination(U::Universe, X::LazySet)
    return _linear_combination_universe(U, X)
end

# ============== #
# disambiguation #
# ============== #

for T in (:ConvexSet, :Universe)
    @eval @commutative function linear_combination(∅::EmptySet, X::($T))
        return _linear_combination_emptyset(∅, X)
    end
end

@commutative function linear_combination(U::Universe, X::ConvexSet)
    return _linear_combination_universe(U, X)
end
