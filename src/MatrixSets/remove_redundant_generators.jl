"""
    remove_redundant_generators!(MZ::MatrixZonotope; tol=1e-9)

Remove redundant generators from a matrix zonotope.

# Input

- `MZ` -- a matrix zonotope
- `tol` -- a specified tolerance

# Output

A new matrix zonotope with fewer generators, or the same matrix zonotope 
if no generator could be removed.

# Algorithm 

The function discards generators whose absolute entry values do not exceed the tolerance`tol`.
"""
function remove_redundant_generators(MZ::MatrixZonotope; tol=1e-9)
    Gs = generators(MZ)
    idx = indexvector(MZ)
    
    @inbounds for i in reverse(eachindex(Gs))
        if maximum(abs, Gs[i]) < tol
            deleteat!(Gs, i)
            deleteat!(idx, i)
        end
    end
    return MatrixZonotope(center(MZ), Gs, idx)
end
