"""
    an_element(hs::HalfSpace)

Return some element of a half-space.

### Input

- `hs` -- half-space

### Output

An element on the defining hyperplane.
"""
function an_element(hs::HalfSpace)
    return _an_element_helper_hyperplane(hs.a, hs.b)
end
