"""
# Extended help

    an_element(hs::HalfSpace)

### Output

This method computes an element on the defining hyperplane and then moves inside
the half-space region for a robust result.
"""
function an_element(hs::HalfSpace)
    return _an_element_halfspace(hs.a, hs.b; direction_inside=true)
end

# see ext/LazySets/LazySetsHalfSpaceExt.jl
function _an_element_halfspace(a, b;
                               nonzero_entry_a::Int=findfirst(!iszero, a),
                               direction_inside::Bool)
    return error()
end
