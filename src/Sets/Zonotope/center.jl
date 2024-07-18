"""
    center(Z::Zonotope)

Return the center of a zonotope.

### Input

- `Z` -- zonotope

### Output

The center of the zonotope.
"""
function center(Z::Zonotope)
    return Z.center
end
