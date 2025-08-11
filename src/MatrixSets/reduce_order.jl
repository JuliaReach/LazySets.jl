"""
    reduce_order(MZ::MatrixZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

Reduce the order of a zonotopic set by overapproximating with a zonotope with
fewer generators.

### Input

- `MZ`     -- matrix zonotope
- `r`      -- desired order
- `method` -- (optional, default: `GIR05()`) the reduction method used

### Output

A new zonotope with fewer generators, if possible.

# Algorithm 

This function first vectorizes the matrix zonotope into a standard zonotope, 
reduces the order of the resulting zonotope, and then converts it back 
to a matrix zonotope of the original dimensions.

# Extended help

    reduce_order(Z::AbstractZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())
"""
function reduce_order(MZ::MatrixZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    @assert r ≥ 1 "cannot reduce below order 1 (got $r)"

    if order(MZ) <= r
        return MZ
    end

    Z = vectorize(MZ)
    Zred = reduce_order(Z, r, method)

    # reshape to matrix zonotope
    dim = size(MZ)
    return matrixize(Zred, dim)
end
