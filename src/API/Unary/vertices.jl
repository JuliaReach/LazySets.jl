@doc """
    vertices(X::LazySet)

Construct an iterator over the vertices of a polytopic set.

### Input

- `X` -- polytopic set

### Output

An iterator over the vertices of `X`.
"""
vertices(::LazySet)
function vertices end
