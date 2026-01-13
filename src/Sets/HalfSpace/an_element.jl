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

@inline function _an_element_halfspace(a::AbstractVector, b;
                                       nonzero_entry_a::Int=findfirst(!iszero, a),
                                       direction_inside::Bool)
    x = _an_element_helper_hyperplane(a, b, nonzero_entry_a)
    if direction_inside
        x .-= a
    else
        x .+= a
    end
    return x
end
