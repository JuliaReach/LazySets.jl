"""
    isempty(cpa::CartesianProductArray)

Check whether a Cartesian product of a finite number of sets is empty.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cpa::CartesianProductArray)
    return any(isempty, array(cpa))
end
