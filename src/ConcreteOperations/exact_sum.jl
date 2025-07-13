# for non-parametric set types, the exact sum coincides with the Minkowski sum
@validate function exact_sum(X::LazySet, Y::LazySet)
    if isparametrictype(typeof(X)) || isparametrictype(typeof(Y))
        throw(ArgumentError("`exact_sum` for set types $(typeof(X)) and " *
                            "$(typeof(Y)) is not implemented"))
    end

    return minkowski_sum(X, Y)
end
