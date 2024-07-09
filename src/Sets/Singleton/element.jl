"""
    element(S::Singleton)

Return the element of a singleton.

### Input

- `S` -- singleton

### Output

The element of the singleton.
"""
function element(S::Singleton)
    return S.element
end
