struct MatrixZonotopeProduct{N, MAT1<:MatrixZonotope{N}, NM, MAT2<:MatrixZonotope{NM}} <: AbstractMatrixZonotope{N}
    A::MAT1
    B::MAT2

    function MatrixZonotopeProduct(A::MAT1, B::MAT2) where {N, MAT1<:MatrixZonotope{N}, NM, MAT2<:MatrixZonotope{NM}}
        @assert size(A.A0, 2) == size(B.A0, 1) "incompatible dimensions"
        return new{N, MAT1, NM, MAT2}(A, B)
    end
end

@commutative *(A::Real, B::MatrixZonotope) = scale(A, B)

function *(A::MatrixZonotope, B::MatrixZonotope)
    return MatrixZonotopeProduct(A, B)
end

Base.size(M::MatrixZonotopeProduct) = (size(M.A, 1), size(M.B, 2))
Base.size(M::MatrixZonotopeProduct, d::Int) = size(M)[d]
