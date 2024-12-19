"""
# Extended help

    an_element(U::Universe{N}) where {N}

### Algorithm

The output is the origin.
"""
function an_element(U::Universe{N}) where {N}
    return zeros(N, dim(U))
end
