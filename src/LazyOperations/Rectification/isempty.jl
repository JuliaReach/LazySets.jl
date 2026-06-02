"""
    isempty(R::Rectification)

Check whether a rectification is empty.

### Input

- `R` -- rectification

### Output

`true` iff the wrapped set is empty.
"""
function isempty(R::Rectification)
    return isempty(R.X)
end
