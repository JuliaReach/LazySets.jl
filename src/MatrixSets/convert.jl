function load_intervalmatrices_conversion()
    return quote
        using .IntervalMatrices: IntervalMatrix, mid, sup, radius

        """
            convert(::Type{MatrixZonotope}, IM::IntervalMatrix)

        Convert an interval matrix to a matrix zonotope

        ### Input

        - `MatrixZonotope` -- target type
        - `IM` -- an interval matrix

        ### Output

        A matrix zonotope with one generator

        ### Example 

        ```jldoctest
        julia> using LazySets, IntervalMatrices

        julia> IM = IntervalMatrix([interval(-1.1, -0.9) interval(-4.1, -3.9);
                    interval(3.9, 4.1) interval(-1.1, -0.9)])
        2×2 IntervalMatrix{Float64, IntervalArithmetic.Interval{Float64}, Matrix{IntervalArithmetic.Interval{Float64}}}:
         [-1.10001, -0.9]  [-4.1, -3.89999]
          [3.89999, 4.1]   [-1.10001, -0.9]

        julia> MZ = convert(MatrixZonotope, IM)
        MatrixZonotope{Float64, Matrix{Float64}}([-1.0 -4.0; 4.0 -1.0], [[0.10000000000000009 0.0; 0.0 0.0], [0.0 0.0; 0.10000000000000009 0.0], [0.0 0.10000000000000009; 0.0 0.0], [0.0 0.0; 0.0 0.10000000000000009]], [1, 2, 3, 4])
        ```
        """
        function Base.convert(::Type{MatrixZonotope}, IM::IntervalMatrix{T}) where T
            #TODO change to sparse matrix -> problem: all the generators of a MZ should have the same type
            n, m = size(IM)
            
            # center must be the midpoint matrix
            center = mid(IM)

            # radius matrix gives the half-widths
            halfIM = radius(IM)

            gens = Vector{Matrix{T}}(undef, n * m)
            idx = 1
            @inbounds for j in 1:m
                for i in 1:n
                    G = zeros(T, n, m)
                    G[i, j] = halfIM[i, j]
                    gens[idx] = G
                    idx += 1
                end
            end

            return MatrixZonotope(center, gens)
        end
    end
end
