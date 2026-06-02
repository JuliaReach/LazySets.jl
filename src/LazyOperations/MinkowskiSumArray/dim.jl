"""
    dim(msa::MinkowskiSumArray)

Return the dimension of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

The ambient dimension of the Minkowski sum of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(msa::MinkowskiSumArray)
    return length(msa.array) == 0 ? 0 : dim(msa.array[1])
end
