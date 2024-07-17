"""
    iscomplement(H1::HalfSpace{N}, H2::HalfSpace) where {N}

Check if two half-spaces complement each other.

### Input

- `H1` -- half-space
- `H2` -- half-space

### Output

`true` iff `H1` and `H2` are complementary, i.e., have opposite normal
directions and identical boundaries (defining hyperplanes).
"""
function iscomplement(H1::HalfSpace{N}, H2::HalfSpace) where {N}
    # check that the half-spaces have converse directions
    res, factor = ismultiple(H1.a, H2.a)
    if !res || !_leq(factor, zero(N))
        return false
    end

    # check that the half-spaces touch each other
    return _isapprox(H1.b, factor * H2.b)
end
