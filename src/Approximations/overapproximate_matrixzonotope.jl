"""
    overapproximate(MZP::MatrixZonotopeProduct{N,S},
                         ::Type{<:MatrixZonotope}) where {N,S<:AbstractMatrix{N}}

Overapproximate the product of matrix zonotopes, following [AlthoffKS11; Equation 4.10](@citet).

### Input

- `MZP` -- a `MatrixZonotopeProduct`
- `MatrixZonotope` -- target type

### Output

A matrix zonotope overapproximating the matrix zonotope product
"""
function overapproximate(MZP::MatrixZonotopeProduct{N,S},
                         ::Type{<:MatrixZonotope}) where {N,S<:AbstractMatrix{N}}
    facs = factors(MZP)
    if length(facs) == 1
        return facs[1]
    end

    return foldl(facs) do A, B
        nB = ngens(B)
        nA = ngens(A)

        if nA == 0
            return linear_map(center(A), B)
        elseif nB == 0
            return linear_map(A, center(B))
        end

        gens = Vector{S}(undef, nB + nA + nA * nB)
        idx = 1

        # G₀ * Hⱼ
        @inbounds for Hj in generators(B)
            gens[idx] = center(A) * Hj
            idx += 1
        end

        # Gᵢ * H₀
        @inbounds for Gi in generators(A)
            gens[idx] = Gi * center(B)
            idx += 1
        end

        # Gᵢ * Hⱼ
        @inbounds for Gi in generators(A), Hj in generators(B)
            gens[idx] = Gi * Hj
            idx += 1
        end

        c = center(A) * center(B)
        res = MatrixZonotope(c, gens)
        return remove_redundant_generators(res)
    end
end
