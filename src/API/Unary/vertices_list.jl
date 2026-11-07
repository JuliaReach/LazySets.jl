@doc """
    vertices_list(X::LazySet)

Compute a list of vertices of a polytopic set.

### Input

- `X` -- polytopic set

### Output

A list of the vertices of `X`.
"""
vertices_list(::LazySet)
function vertices_list end
