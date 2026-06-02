"""
    ρ(d::AbstractVector, msa::MinkowskiSumArray)

Evaluate the support function of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum of a finite number of sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the Minkowski sum of multiple sets evaluations to the
sum of the support-function evaluations of each set.
"""
@validate function ρ(d::AbstractVector, msa::MinkowskiSumArray)
    return sum(ρ(d, Xi) for Xi in msa.array)
end
