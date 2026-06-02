"""
    an_element(tr::Translation)

Return some element of a translation.

### Input

- `tr` -- translation of a set

### Output

An element in the translation.

### Notes

This function first asks for `an_element` of the wrapped set, then translates
this element according to the given translation vector.
"""
function an_element(tr::Translation)
    return an_element(tr.X) + tr.v
end
