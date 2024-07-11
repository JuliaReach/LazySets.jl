"""
    direction(L::Line)

Return the direction of the line.

### Input

- `L` -- line

### Output

The direction of the line.

### Notes

The direction is not necessarily normalized.
See [`normalize(::Line, ::Real)`](@ref) / [`normalize!(::Line, ::Real)`](@ref).
"""
direction(L::Line) = L.d
