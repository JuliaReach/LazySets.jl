"""
    isbounded(cpa::CartesianProductArray)

Check whether a Cartesian product of a finite number of sets is bounded.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cpa::CartesianProductArray)
    return all(isbounded, cpa.array)
end
