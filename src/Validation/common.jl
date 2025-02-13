function validate_dim(X::LazySet, n::Int; fun::Function)
    if dim(X) != n
        throw(DimensionMismatch("`$(string(fun))` requires a $n-dimensional " *
                                "set but received a $(dim(X))-dimensional set"))
    end
    return true
end

function validate_dims(X::LazySet, dims; fun::Function)
    if dim(X) ∉ dims
        throw(DimensionMismatch("`$(string(fun))` requires a set of dimensions " *
                                "$dims but received a $(dim(X))-dimensional set"))
    end
    return true
end

function validate_same_dim(X::LazySet, Y::LazySet; fun::Function)
    if dim(X) != dim(Y)
        throw(DimensionMismatch("`$(string(fun))` requires sets of the same dimension but " *
                                "received sets of dimensions $(dim(X)) and $(dim(Y))"))
    end
    return true
end
