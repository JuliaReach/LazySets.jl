"""
    ρ(d::AbstractVector, cpa::CartesianProductArray)

Evaluate the support function of a Cartesian product of a finite number of sets.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product of a finite number of sets

### Output

The evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
@validate function ρ(d::AbstractVector, cpa::CartesianProductArray)
    N = promote_type(eltype(d), eltype(cpa))
    sfun = zero(N)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        sfun += ρ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return sfun
end

# faster version for sparse vectors
@validate function ρ(d::AbstractSparseVector, cpa::CartesianProductArray)
    N = promote_type(eltype(d), eltype(cpa))
    # idea: see the σ method for AbstractSparseVector
    sfun = zero(N)
    indices, _ = findnz(d)
    if isempty(indices)
        # direction is the zero vector
        return sfun
    end
    next_idx = 1
    next_dim = indices[next_idx]
    m = length(indices)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if next_dim <= i1
            # there is a non-zero entry in this block
            sfun += ρ(d[i0:i1], Xi)

            # find next index outside the current block
            next_idx += 1
            while next_idx <= m && indices[next_idx] <= i1
                next_idx += 1
            end
            if next_idx > m
                # no more non-zero entries
                break
            end
            next_dim = indices[next_idx]
        end
        i0 = i1 + 1
    end
    return sfun
end

# faster version for single-entry vectors
@validate function ρ(d::SingleEntryVector, cpa::CartesianProductArray)
    i0 = 1
    idx = d.i
    for Xi in cpa.array
        ni = dim(Xi)
        i1 = i0 + ni - 1
        if i0 <= idx && idx <= i1
            return ρ(SingleEntryVector(d.i - i0 + 1, ni, d.v), Xi)
        end
        i0 = i1 + 1
    end
    return sfun
end
