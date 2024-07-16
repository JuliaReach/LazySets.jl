"""
    ∈(v::AbstractVector, X::Star)

Check whether a given point is contained in a star.

### Input

- `v` -- point/vector
- `X` -- star

### Output

`true` iff ``v ∈ X``.

### Algorithm

The implementation is identical to
[`∈(::AbstractVector, ::AbstractAffineMap)`](@ref).
"""
function ∈(x::AbstractVector, X::Star)
    return basis(X) \ (x - center(X)) ∈ predicate(X)
end
