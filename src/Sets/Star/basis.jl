"""
    basis(X::Star)

Return the basis vectors of a star.

### Input

- `X` -- star

### Output

A matrix where each column is a basis vector of the star.
"""
basis(X::Star) = X.V
