function load_intervalmatrices_conversion()
    return quote
        using .IntervalMatrices: IntervalMatrix, mid, sup

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
        MatrixZonotope{Float64, Matrix{Float64}}([-1.0 -4.0; 4.0 -1.0],
        [[0.09999999999999998 0.10000000000000009; 0.09999999999999964 0.09999999999999998]], [1])
        ````
        """
        function Base.convert(::Type{MatrixZonotope}, IM::IntervalMatrix)
            c = mid(IM)
            G = [radius(IM)]
            return MatrixZonotope(c, G)
        end
    end
end
