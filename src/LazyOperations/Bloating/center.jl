"""
    center(B::Bloating)

Return the center of a bloated set.

### Input

- `B` -- bloated set

### Output

The center of the wrapped set.

### Notes

This implementation disregards negative bloating, which could potentially remove
the center from the set.
"""
function center(B::Bloating)
    return center(B.X)
end
