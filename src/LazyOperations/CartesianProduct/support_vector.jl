"""
    σ(d::AbstractVector, cp::CartesianProduct)

Return a support vector of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
@validate function σ(d::AbstractVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    return [σ(d[1:n1], cp.X); σ(d[(n1 + 1):length(d)], cp.Y)]
end

# faster version for single-entry vectors
@validate function σ(d::SingleEntryVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    idx = d.i
    if idx <= n1
        return [σ(SingleEntryVector(idx, n1, d.v), cp.X); an_element(cp.Y)]
    else
        n2 = length(d) - n1
        return [an_element(cp.X); σ(SingleEntryVector(idx - n1, n2, d.v), cp.Y)]
    end
end
