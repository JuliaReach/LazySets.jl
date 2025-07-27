"""
    MatrixZonotopeProduct{N}(factors::Vector{<:MatrixZonotope{N}}) <: AbstractMatrixZonotope{N}

Represents the product of multiple matrix zonotopes.

### Fields

- `factors` -- a vector of matrix zonotopes `[A1, A2, ..., An]`

### Notes

Mathematically, this represents the set:

```math
\\mathcal{C} = \\{ A_1 A_2 \\dots A_n ~|~ A_i \\in \\mathcal{A}_i \\}
```
"""
struct MatrixZonotopeProduct{N,S<:AbstractMatrix{N}} <: AbstractMatrixZonotope{N}
    factors::Vector{MatrixZonotope{N,S}}

    function MatrixZonotopeProduct(factors::Vector{<:MatrixZonotope{N,T}}) where {N,
                                                                                  T<:AbstractMatrix{N}}
        for i in 1:(length(factors) - 1)
            @assert size(factors[i], 2) == size(factors[i + 1], 1) "incompatible dimensions for factors at index $i ($(size(factors[i],2))) and $(i+1) ($(size(factors[i+1],1)))"
        end
        return new{N,T}(factors)
    end
end

MatrixZonotopeProduct(ms::MatrixZonotope...) = MatrixZonotopeProduct(collect(ms))

Base.:*(A::MatrixZonotope, B::MatrixZonotope) = MatrixZonotopeProduct([A, B])
Base.:*(P::MatrixZonotopeProduct, B::MatrixZonotope) = MatrixZonotopeProduct(vcat(P.factors, B))
Base.:*(A::MatrixZonotope, P::MatrixZonotopeProduct) = MatrixZonotopeProduct(vcat(A, P.factors))
Base.:*(P1::MatrixZonotopeProduct, P2::MatrixZonotopeProduct) = MatrixZonotopeProduct(vcat(P1.factors, P2.factors))

"""
    factors(MZP::MatrixZonotopeProduct)

Return the factors of a matrix zonotope product.
"""
function factors(MZP::MatrixZonotopeProduct)
    return MZP.factors
end

"""
    nfactors(MZP::MatrixZonotopeProduct)

Return the number of factors of a matrix zonotope product.
"""
function nfactors(MZP::MatrixZonotopeProduct)
    return length(factors(MZP))
end

Base.size(P::MatrixZonotopeProduct) = (size(P.factors[1].A0, 1), size(P.factors[end].A0, 2))
Base.size(P::MatrixZonotopeProduct, d::Int) = size(P)[d]

Base.:(==)(P1::MatrixZonotopeProduct, P2::MatrixZonotopeProduct)= P1.factors == P2.factors
