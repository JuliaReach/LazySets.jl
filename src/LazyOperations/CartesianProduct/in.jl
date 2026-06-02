"""
    in(x::AbstractVector, cp::CartesianProduct)

Check whether a given point is contained in a Cartesian product.

### Input

- `x`  -- point/vector
- `cp` -- Cartesian product

### Output

`true` iff ``x ∈ cp``.
"""
@validate function in(x::AbstractVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    return view(x, 1:n1) ∈ cp.X &&
           view(x, (n1 + 1):length(x)) ∈ cp.Y
end
