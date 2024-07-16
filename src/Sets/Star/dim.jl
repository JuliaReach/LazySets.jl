"""
    dim(X::Star)

Return the dimension of a star.

### Input

- `X` -- star

### Output

The ambient dimension of a star.
"""
function dim(X::Star)
    return length(X.c)
end
