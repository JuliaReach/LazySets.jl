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

# disambiguation
@commutative function linear_combination(∅::EmptySet, X::ConvexSet)
    return _linear_combination_emptyset(∅, X)
end

function linear_combination(P1::AbstractPolynomialZonotope,
                            P2::AbstractPolynomialZonotope)
    SSPZ1 = convert(SimpleSparsePolynomialZonotope, P1)
    SSPZ2 = convert(SimpleSparsePolynomialZonotope, P2)
    return linear_combination(SSPZ1, SSPZ2)
end

@commutative function linear_combination(Z::AbstractZonotope, P::AbstractPolynomialZonotope)
    return linear_combination(convert(SparsePolynomialZonotope, Z), P)
end
