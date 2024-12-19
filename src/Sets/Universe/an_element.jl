"""
# Extended help

    an_element(U::Universe)

### Algorithm

The output is the origin.
"""
function an_element(U::Universe)
    N = eltype(U)
    return zeros(N, dim(U))
end
