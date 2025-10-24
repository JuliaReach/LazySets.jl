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

function validate_map_dim(M::AbstractMatrix, X::LazySet; fun::Function)
    if size(M, 2) != dim(X)
        throw(DimensionMismatch("`$(string(fun))` requires a matrix and a set of compatible " *
                                "dimension but received a matrix of size $(size(M)) and a " *
                                "set of dimension $(dim(X))"))
    end
    return true
end

function validate_index(i::Int, X::LazySet; fun::Function)
    if i < 1 || i > dim(X)
        throw(DimensionMismatch("`$(string(fun))` requires an index for " *
                                "dimension $(dim(X)) but received $i"))
    end
    return true
end

function validate_index_vector(v::AbstractVector{Int}, X::LazySet; fun::Function)
    n = dim(X)
    b = falses(n)
    for e in v
        if e < 1 || e > n
            throw(DimensionMismatch("`$(string(fun))` requires an index vector " *
                                    "for dimension $n but received the vector $v"))
        end
        if b[e]
            throw(ArgumentError("an index vector must not contain duplicate entries"))
        else
            b[e] = true
        end
    end
    return true
end

function validate_pnorm(p::Real; fun::Function)
    if p < one(p)
        throw(ArgumentError("`$(string(fun))` requires a p-norm for p ≥ 1 " *
                            "but received $(p)"))
    end
    return true
end
