"""
    normalize(L::Line{N}, p::Real=N(2)) where {N}

Normalize the direction of a line.

### Input

- `L` -- line
- `p` -- (optional, default: `2.0`) vector `p`-norm used in the normalization

### Output

A line whose direction has unit norm w.r.t. the given `p`-norm.

### Notes

See also [`normalize!(::Line, ::Real)`](@ref) for the in-place version.
"""
function normalize(L::Line{N}, p::Real=N(2)) where {N}
    return normalize!(copy(L), p)
end

"""
    normalize!(L::Line{N}, p::Real=N(2)) where {N}

Normalize the direction of a line storing the result in `L`.

### Input

- `L` -- line
- `p` -- (optional, default: `2.0`) vector `p`-norm used in the normalization

### Output

A line whose direction has unit norm w.r.t. the given `p`-norm.
"""
function normalize!(L::Line{N}, p::Real=N(2)) where {N}
    normalize!(L.d, p)
    return L
end
