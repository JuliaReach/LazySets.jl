"""
    center(cpa::CartesianProductArray)

Compute the center of a Cartesian product of a finite number of
centrally-symmetric sets.

### Input

- `cpa` -- Cartesian product of a finite number of centrally-symmetric sets

### Output

The center of the Cartesian product of a finite number of sets.
"""
function center(cpa::CartesianProductArray)
    return reduce(vcat, center(X) for X in cpa)
end
