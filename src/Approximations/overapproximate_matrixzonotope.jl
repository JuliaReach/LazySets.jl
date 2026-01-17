function load_intervalmatrices_overapproximation_matrixzonotope()
    return quote
        using .IntervalMatrices: IntervalMatrix

        """
            overapproximate(expA::MatrixZonotopeExp{N,T}, ::Type{<:MatrixZonotope},
                                 k::Int=2) where {N,T<:AbstractMatrixZonotope{N}}

        Overapproximate the matrix zonotope exponential ``exp(\\mathcal{A})``

        ### Input

        - `expA` -- `MatrixZonotopeExp`
        - `MatrixZonotope` -- target type
        - `k` -- (default: `2`) the order of the Taylor expansion
        - `tol` -- (default: `1e-9`) tolerance used when pruning generators after the overapproximation

        ### Output

        A matrix zonotope overapproximating the matrix zonotope exponential

        ### Algorithm

        The expansion

        ```math
        exp(\\mathcal{A}) ⊆ \\sum_i^k \\frac{\\mathcal{A}^i}{i!} + E_k
        ```

        is computed by overapproximating the matrix zonotope powers ``A^i``
        for ``i=0, …, k``.
        The remainder term ``E_k`` is computed through interval arithmetic
        following [AlthoffKS11; Proposition 4.1](@citet).
        """
        function overapproximate(expA::MatrixZonotopeExp{N,T},
                                 ::Type{<:MatrixZonotope},
                                 k::Int=2; tol::Real=1e-9) where {N,T<:AbstractMatrixZonotope{N}}

            # overapproximate the product MZP = A*B*... ---
            MZP = MatrixZonotopeProduct(expA.M)
            MZ = overapproximate(MZP, MatrixZonotope)

            Id = MatrixZonotope(Matrix{N}(I, size(MZ)), Matrix{N}[])
            powers = Vector{typeof(MZ)}(undef, k)
            powers[1] = MZ
            @inbounds for i in 2:k
                term = overapproximate(MZ * powers[i - 1], MatrixZonotope)
                powers[i] = scale(1 / i, term)
            end
            W = reduce(minkowski_sum, powers)
            W = minkowski_sum(W, Id)

            IM = overapproximate(MZ, IntervalMatrix)
            E = IntervalMatrices._exp_remainder(IM, N(1), k)
            res = minkowski_sum(W, E)

            return remove_redundant_generators(res; tol=tol)
        end

        # TODO avoid `convert`
        function LazySets.minkowski_sum(MZ::MatrixZonotope, IM::IntervalMatrix)
            return minkowski_sum(MZ, convert(MatrixZonotope, IM))
        end

        """
            overapproximate(A::MatrixZonotope, ::Type{<:IntervalMatrix})

        Overapproximate a matrix zonotope with an interval matrix.

        ### Input

        - `A` -- a matrix zonotope
        - `IntervalMatrix` -- target type

        ### Output

        An interval matrix overapproximating the matrix zonotope.
        """
        function overapproximate(A::MatrixZonotope, ::Type{<:IntervalMatrix})
            A_abs = sum([abs.(Ai) for Ai in generators(A)])
            A₊ = center(A) + A_abs
            A₋ = center(A) - A_abs
            return IntervalMatrix(A₋, A₊)
        end
    end
end


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
