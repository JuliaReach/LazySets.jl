"""
# Extended help

    isuniversal(hs::HalfSpace, [witness]::Bool=false)

### Algorithm

Witness production works analogously to `an_element`.
"""
function isuniversal(hs::HalfSpace, witness::Bool=false)
    if witness
        v = _non_element_halfspace(hs.a, hs.b)
        return (false, v)
    else
        return false
    end
end

function _non_element_halfspace(a, b)
    return _an_element_halfspace(a, b; direction_inside=false)
end
