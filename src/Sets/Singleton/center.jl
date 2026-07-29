"""
    center(S::Singleton)

Return the center of a singleton.

### Input

- `S` -- singleton

### Output

The unique element of the singleton.
"""
function center(S::Singleton)
    return S.element
end
