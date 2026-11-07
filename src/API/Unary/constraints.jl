@doc """
    constraints(X::LazySet)

Construct an iterator over the constraints of a polyhedral set.

### Input

- `X` -- polyhedral set

### Output

An iterator over the constraints of `X`.
"""
constraints(::LazySet)
function constraints end
