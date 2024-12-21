"""
# Extended help

    isuniversal(H::Hyperplane, [witness]::Bool=false)

### Algorithm

A witness is produced by adding the normal vector to an element on the
hyperplane.
"""
function isuniversal(H::Hyperplane, witness::Bool=false)
    if witness
        v = _non_element_halfspace(H.a, H.b)
        return (false, v)
    else
        return false
    end
end
