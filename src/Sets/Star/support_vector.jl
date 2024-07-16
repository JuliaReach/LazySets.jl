"""
    σ(d::AbstractVector, X::Star)

Return a support vector of a star.

### Input

- `d` -- direction
- `X` -- star

### Output

A support vector in the given direction.
"""
function σ(d::AbstractVector, X::Star)
    A = basis(X)
    return A * σ(At_mul_B(A, d), predicate(X)) + center(X)
end
