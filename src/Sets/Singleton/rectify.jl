"""
    rectify(S::Singleton)

Concrete rectification of a singleton.

### Input

- `S` -- singleton

### Output

The `Singleton` that corresponds to the rectification of `S`.
"""
function rectify(S::Singleton)
    return Singleton(rectify(element(S)))
end
