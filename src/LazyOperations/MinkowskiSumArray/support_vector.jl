"""
    σ(d::AbstractVector, msa::MinkowskiSumArray)

Return a support vector of a Minkowski sum of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
@validate function σ(d::AbstractVector, msa::MinkowskiSumArray)
    return _σ_msum_array(d, msa.array)
end

@inline function _σ_msum_array(d::AbstractVector, array::AbstractVector{<:LazySet})
    return sum(σ(d, Xi) for Xi in array)
end
