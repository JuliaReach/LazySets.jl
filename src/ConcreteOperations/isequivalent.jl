"""
# Extended help

    isequivalent(X::LazySet, Y::LazySet)

### Algorithm

The default implementation first checks `X ≈ Y`, which returns `true` if and
only if `X` and `Y` have the same base type and approximately the same values.
If that fails, the implementation checks the double inclusion `X ⊆ Y && Y ⊆ X`.

### Examples

```jldoctest
julia> X = BallInf([0.1, 0.2], 0.3);

julia> Y = convert(HPolytope, X);

julia> X == Y
false

julia> isequivalent(X, Y)
true
```
"""
function isequivalent(X::LazySet, Y::LazySet)
    try  # TODO temporary try-catch construct until ≈ is fixed for all set types
        if X ≈ Y
            return true
        end
    catch e
    end
    return _isequivalent_inclusion(X, Y)
end

function _isequivalent_inclusion(X::LazySet, Y::LazySet)
    return X ⊆ Y && Y ⊆ X
end

@commutative function isequivalent(S::AbstractSingleton, Z::AbstractZonotope)
    return center(Z) == element(S) && iszero(genmat(Z))
end

function isequivalent(S1::AbstractSingleton, S2::AbstractSingleton)
    return isapprox(element(S1), element(S2))
end
