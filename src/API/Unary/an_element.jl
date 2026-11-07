@doc """
    an_element(X::LazySet)

Return some element of a nonempty set.

### Input

- `X` -- set

### Output

An element of `X` unless `X` is empty.
"""
an_element(::LazySet)
function an_element end
