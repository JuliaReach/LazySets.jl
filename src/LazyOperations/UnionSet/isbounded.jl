"""
    isbounded(cup::UnionSet)

Check whether the union of two sets is bounded.

### Input

- `cup` -- union of two sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSet)
    return isbounded(cup.X) && isbounded(cup.Y)
end
