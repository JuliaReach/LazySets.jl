# the only parametric sets in this library are polynomial zonotopes;
# hence, we compute the normal Minkowski sum for all other set combinations
function exact_sum(X::LazySet, Y::LazySet)
    @assert dim(X) == dim(Y) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(Y)), respectively"

    if X isa AbstractPolynomialZonotope || Y isa AbstractPolynomialZonotope
        throw(ArgumentError("`exact_sum` for set types $(typeof(X)) and " *
                            "$(typeof(Y)) is not implemented"))
    end
    return minkowski_sum(X, Y)
end
