"""
    center(sih::SymmetricIntervalHull, i::Int)

Return the center along a given dimension of the symmetric interval hull of a
set.

### Input

- `sih` -- symmetric interval hull of a set
- `i`   -- dimension of interest

### Output

The center along a given dimension of the symmetric interval hull of a set.
"""
@validate function center(sih::SymmetricIntervalHull, i::Int)
    N = eltype(sih)
    return zero(N)
end

"""
    center(sih::SymmetricIntervalHull)

Return the center of the symmetric interval hull of a set.

### Input

- `sih` -- symmetric interval hull of a set

### Output

The origin.
"""
function center(sih::SymmetricIntervalHull)
    N = eltype(sih)
    return zeros(N, dim(sih))
end
