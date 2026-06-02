"""
    isbounded(cup::UnionSetArray)

Check whether the union of a finite number of sets is bounded.

### Input

- `cup` -- union of a finite number of sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSetArray)
    return all(isbounded, array(cup))
end
