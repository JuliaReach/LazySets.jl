"""
    dim(sih::SymmetricIntervalHull)

Return the dimension of the symmetric interval hull of a set.

### Input

- `sih` -- symmetric interval hull of a set

### Output

The ambient dimension of the symmetric interval hull of a set.
"""
function dim(sih::SymmetricIntervalHull)
    return dim(sih.X)
end
