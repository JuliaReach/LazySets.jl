"""
# Extended help

    ⊂(X::LazySet, Y::LazySet, [witness]::Bool=false)

### Algorithm

The default implementation first checks inclusion of `X` in `Y` and then checks
noninclusion of `Y` in `X`:
"""
@validate function ⊂(X::LazySet, Y::LazySet, witness::Bool=false)
    @assert dim(X) == dim(Y) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(Y)), respectively"

    if witness
        N = promote_type(eltype(X), eltype(Y))
        res, w = ⊆(X, Y, witness)
        if res
            res, w = ⊆(Y, X, witness)
            if res
                return (false, N[])
            else
                return (true, w)
            end
        else
            return (false, N[])
        end
    end
    return (X ⊆ Y) && !(Y ⊆ X)
end
