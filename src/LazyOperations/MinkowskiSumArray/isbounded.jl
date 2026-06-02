"""
    isbounded(msa::MinkowskiSumArray)

Check whether a Minkowski sum of a finite number of sets is bounded.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(msa::MinkowskiSumArray)
    return all(isbounded, msa.array)
end
