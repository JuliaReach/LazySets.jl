"""
# Extended help

    isuniversal(L::Line2D, [witness]::Bool=false)

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(L::Line2D, witness::Bool=false)
    if witness
        v = _non_element_halfspace(L.a, L.b)
        return (false, v)
    else
        return false
    end
end
