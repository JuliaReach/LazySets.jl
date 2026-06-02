"""
    an_element(R::Rectification)

Return some element of a rectification.

### Input

- `R` -- rectification

### Output

An element in the rectification.
The implementation relies on the `an_element` function of the wrapped set.
"""
function an_element(R::Rectification)
    return rectify(an_element(R.X))
end
