"""
    an_element(∅::EmptySet)

Return some element of an empty set.

### Input

- `∅` -- empty set

### Output

An error.
"""
function an_element(::EmptySet)
    return error("an empty set does not contain any element")
end
