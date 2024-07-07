"""
    element(Z::ZeroSet{N}) where {N}

Return the element of a zero set.

### Input

- `Z` -- zero set

### Output

The element of the zero set, i.e., a zero vector.
"""
function element(Z::ZeroSet{N}) where {N}
    return zeros(N, Z.dim)
end

"""
    element(Z::ZeroSet{N}, ::Int) where {N}

Return the i-th entry of the element of a zero set.

### Input

- `Z` -- zero set
- `i` -- dimension

### Output

The i-th entry of the element of the zero set, i.e., 0.
"""
function element(Z::ZeroSet{N}, ::Int) where {N}
    return zero(N)
end
