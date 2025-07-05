function validate_same_dim(X::LazySet, Y::LazySet; fun::Function)
    if dim(X) != dim(Y)
        throw(DimensionMismatch("`$(string(fun))` requires sets of the same dimension but " *
                                "received sets of dimensions $(dim(X)) and $(dim(Y))"))
    end
    return true
end

function validate_same_dim(x::AbstractVector, X::LazySet; fun::Function)
    if length(x) != dim(X)
        throw(DimensionMismatch("`$(string(fun))` requires a vector and a set of the same " *
                                "dimension but received a vector of length $(length(x)) and a " *
                                "set of dimension $(dim(X))"))
    end
    return true
end

function validate_same_dim(n::Int, X::LazySet; fun::Function)
    if n != dim(X)
        throw(DimensionMismatch("`$(string(fun))` requires a $n-dimensional " *
                                "set but received a $(dim(X))-dimensional set"))
    end
    return true
end

function validate_map_dim(M::AbstractMatrix, X::LazySet; fun::Function)
    if size(M, 2) != dim(X)
        throw(DimensionMismatch("`$(string(fun))` requires a matrix and a set of compatible " *
                                "dimension but received a matrix of size $(size(M)) and a " *
                                "set of dimension $(dim(X))"))
    end
    return true
end

function validate_index_vector(v::AbstractVector{Int}, X::LazySet; fun::Function)
    n = dim(X)
    for e in v
        if e < 1 || e > n
            throw(DimensionMismatch("`$(string(fun))` requires an index vector " *
                                    "for dimension $n but received the vector $v"))
        end
    end
    return true
end
