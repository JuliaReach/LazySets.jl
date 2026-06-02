"""
    dim(R::Rectification)

Return the dimension of a rectification.

### Input

- `R` -- rectification

### Output

The ambient dimension of the rectification.
"""
function dim(R::Rectification)
    return dim(R.X)
end
