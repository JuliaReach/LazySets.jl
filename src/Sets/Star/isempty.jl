"""
    isempty(X::Star)

Check whether a star is empty.

### Input

- `X` -- star

### Output

`true` iff the predicate is empty.
"""
function isempty(X::Star)
    return isempty(predicate(X))
end
