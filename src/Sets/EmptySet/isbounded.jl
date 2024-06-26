"""
    isbounded(∅::EmptySet)

Check whether an empty set is bounded.

### Input

- `∅` -- empty set

### Output

`true`.
"""
function isbounded(::EmptySet)
    return true
end
