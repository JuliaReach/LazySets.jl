"""
    isempty(cup::UnionSet)

Check whether the union of two sets is empty.

### Input

- `cup` -- union of two sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSet)
    return isempty(cup.X) && isempty(cup.Y)
end
