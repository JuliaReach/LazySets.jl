"""
    center(Z::ZeroSet)

Return the center of a zero set.

### Input

- `Z` -- zero set

### Output

The unique element of the zero set, i.e., a zero vector.
"""
function center(Z::ZeroSet)
    N = eltype(Z)
    return zeros(N, Z.dim)
end

"""
    center(Z::ZeroSet, i::Int)

Return the i-th entry of the center of a zero set.

### Input

- `Z` -- zero set
- `i` -- dimension

### Output

The i-th entry of the center of the zero set, i.e., 0.
"""
@validate function center(Z::ZeroSet, i::Int)
    N = eltype(Z)
    return zero(N)
end
