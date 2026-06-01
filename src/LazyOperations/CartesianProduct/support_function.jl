"""
    ρ(d::AbstractVector, cp::CartesianProduct)

Evaluate the support function of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
@validate function ρ(d::AbstractVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    return ρ(d[1:n1], cp.X) + ρ(d[(n1 + 1):length(d)], cp.Y)
end

# faster version for single-entry vectors
@validate function ρ(d::SingleEntryVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    idx = d.i
    if idx <= n1
        return ρ(SingleEntryVector(idx, n1, d.v), cp.X)
    else
        n2 = length(d) - n1
        return ρ(SingleEntryVector(idx - n1, n2, d.v), cp.Y)
    end
end
