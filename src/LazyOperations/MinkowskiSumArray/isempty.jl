"""
    isempty(msa::MinkowskiSumArray)

Check whether a Minkowski sum of a finite number of sets is empty.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

`true` iff any of the wrapped sets is empty.
"""
function isempty(msa::MinkowskiSumArray)
    return any(isempty, array(msa))
end
