"""
    volume(cpa::CartesianProductArray)

Compute the volume of a Cartesian product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

The volume.
"""
function volume(cpa::CartesianProductArray)
    return prod(volume, array(cpa))
end
