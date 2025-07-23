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
        @assert length(factors) >= 2 "need at least two matrix zonotope factors"
        for i in 1:(length(factors) - 1)
            @assert size(factors[i], 2) == size(factors[i + 1], 1) "incompatible dimensions for factors at index $i ($(size(factors[i],2))) and $(i+1) ($(size(factors[i+1],1)))"
        end
        return new{N,T}(factors)
    end
end

MatrixZonotopeProduct(ms::MatrixZonotope...) = MatrixZonotopeProduct(collect(ms))

"""
    *(A::MatrixZonotope{N,S}, B::MatrixZonotope{N,S}) where {N,S}

Alias to create a `MatrixZonotopeProduct` object.
"""
Base.:*(A::MatrixZonotope, B::MatrixZonotope) = MatrixZonotopeProduct([A, B])
Base.:*(P::MatrixZonotopeProduct, B::MatrixZonotope) = MatrixZonotopeProduct(vcat(P.factors, B))
Base.:*(A::MatrixZonotope, P::MatrixZonotopeProduct) =MatrixZonotopeProduct(vcat(A, P.factors))
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

"""
    remove_redundant_factors(MZP::MatrixZonotopeProduct)

Returns a simplified `MatrixZonotopeProduct` where generator-free factors are absorbed
into adjacent factors via linear maps.

### Input

- `MZP` -- matrix zonotope product.

### Output

A matrix zonotope product with redundant constant factors removed.
"""
function remove_redundant_factors(MZP::MatrixZonotopeProduct)
    factors_ = factors(MZP)
    gens = MatrixZonotope{eltype(factors_[1]), typeof(factors_[1].A0)}[]

    i = 1
    while i < length(factors_)
        MZ = factors_[i]
        next = factors_[i + 1]

        if isempty(generators(MZ))
            push!(gens, linear_map(center(MZ), next))
            i += 2  # skip next since it's been merged
        else
            push!(gens, MZ)
            i += 1
        end
    end

    # in the the last element apply linear map to the left
    if i == length(factors_)
        last = factors_[end]
        if isempty(generators(last)) && !isempty(gens)
            gens[end] = linear_map(gens[end], center(last))
        else
            push!(gens, last)
        end
    end

    return MatrixZonotopeProduct(gens)
end
