"""
    σ(d::AbstractVector, cpa::CartesianProductArray)

Compute a support vector of a Cartesian product of a finite number of sets.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.
"""
@validate function σ(d::AbstractVector, cpa::CartesianProductArray)
    svec = similar(d)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        svec[i0:i1] = σ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return svec
end

# faster version for sparse vectors
@validate function σ(d::AbstractSparseVector, cpa::CartesianProductArray)
    # idea: We walk through the blocks of `cpa` (i.e., the sets `Xi`) and search
    # for corresponding non-zero entries in `d` (stored in `indices`).
    # `next_idx` is the next index of `indices` such that
    # `next_dim = indices[next_idx]` lies in the next block to consider
    # (potentially skipping some blocks).
    svec = similar(d)
    indices, _ = findnz(d)
    if isempty(indices)
        # direction is the zero vector
        return an_element(cpa)
    end
    next_idx = 1
    next_dim = indices[next_idx]
    m = length(indices)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if next_dim <= i1
            # there is a non-zero entry in this block
            svec[i0:i1] = σ(d[i0:i1], Xi)

            # find next index outside the current block
            next_idx += 1
            while next_idx <= m && indices[next_idx] <= i1
                next_idx += 1
            end
            if next_idx <= m
                next_dim = indices[next_idx]
            end
        else
            svec[i0:i1] = an_element(Xi)
        end
        i0 = i1 + 1
    end
    return svec
end

# faster version for single-entry vectors
@validate function σ(d::SingleEntryVector, cpa::CartesianProductArray)
    svec = similar(d)
    i0 = 1
    idx = d.i
    for Xi in cpa.array
        ni = dim(Xi)
        i1 = i0 + ni - 1
        if i0 <= idx && idx <= i1
            svec[i0:i1] = σ(SingleEntryVector(d.i - i0 + 1, ni, d.v), Xi)
        else
            svec[i0:i1] = an_element(Xi)
        end
        i0 = i1 + 1
    end
    return svec
end
