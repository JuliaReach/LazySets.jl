"""
    predicate(X::Star)

Return the predicate of a star.

### Input

- `X` -- star

### Output

A polyhedral set representing the predicate of the star.
"""
predicate(X::Star) = X.P
