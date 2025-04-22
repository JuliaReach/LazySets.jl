"""
   genmat(Z::ZonotopeMD)

Return the generator matrix of a ZonotopeMD.

### Input

- `Z` -- zonotope in normal form

### Output

A matrix where each column represents one generator of the zonotope `Z`.
"""
function genmat(Z::ZonotopeMD)
    D = spdiagm(0 => Z.d)
    return hcat(Z.M, D)
end
