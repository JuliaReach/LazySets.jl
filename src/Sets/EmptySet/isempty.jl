"""
    isempty(∅::EmptySet)

Check if the empty set is empty.

### Input

- `∅` -- empty set

### Output

`true`.
"""
function isempty(::EmptySet)
    return true
end
