"""
    center(msa::MinkowskiSumArray)

Return the center of a Minkowski sum of a finite number of centrally-symmetric
sets.

### Input

- `msa` -- Minkowski sum of a finite number of centrally-symmetric sets

### Output

The center of the set.
"""
function center(msa::MinkowskiSumArray)
    return sum(center(X) for X in msa)
end
