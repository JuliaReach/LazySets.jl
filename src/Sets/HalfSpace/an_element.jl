"""
# Extended help

    an_element(hs::HalfSpace)

### Output

This method returns an element on the defining hyperplane.
"""
function an_element(hs::HalfSpace)
    return _an_element_helper_hyperplane(hs.a, hs.b)
end
