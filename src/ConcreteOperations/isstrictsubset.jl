"""
# Extended help

    ⊂(X::LazySet, Y::LazySet, [witness]::Bool=false)

### Algorithm

The default implementation first checks inclusion of `X` in `Y` and then checks
inclusion of `Y` in `X`:
"""
@validate function ⊂(X::LazySet, Y::LazySet, witness::Bool=false)
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
            return (false, w)
        end
    end
    return (X ⊆ Y) && !(Y ⊆ X)
end

@validate function ⊂(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _isstrictsubset_emptyset(∅, X, witness)
end

@validate function ⊂(X::LazySet, ∅::EmptySet, witness::Bool=false)
    return _isstrictsubset_emptyset2(X, ∅, witness)
end
