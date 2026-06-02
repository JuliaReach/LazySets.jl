"""
    dim(cpa::CartesianProductArray)

Return the dimension of a Cartesian product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

The ambient dimension of the Cartesian product of a finite number of sets, or
`0` if there is no set in the array.
"""
function dim(cpa::CartesianProductArray)
    return length(cpa.array) == 0 ? 0 : sum(dim(Xi) for Xi in cpa.array)
end
