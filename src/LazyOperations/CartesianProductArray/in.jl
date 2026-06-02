"""
    in(x::AbstractVector, cpa::CartesianProductArray)

Check whether a given point is contained in a Cartesian product of a finite
number of sets.

### Input

- `x`   -- point/vector
- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff ``x ∈ \\text{cpa}``.
"""
@validate function in(x::AbstractVector, cpa::CartesianProductArray)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if x[i0:i1] ∉ Xi
            return false
        end
        i0 = i1 + 1
    end
    return true
end
