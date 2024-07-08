"""
    singleton_list(S::Singleton)

Return the vertices of a singleton as a list of singletons.

### Input

- `S` -- singleton

### Output

The list of vertices of `S`, as `Singleton`.
"""
function singleton_list(S::Singleton)
    return [S]
end
