"""
    MatrixZonotopeProduct{N}(factors::Vector{<:MatrixZonotope{N}}) <: AbstractMatrixZonotope{N}

Represents the product of multiple matrix zonotopes.

### Fields

- `factors` -- a vector of matrix zonotopes `[A1, A2, ..., An]`

### Notes

Mathematically, this represents the set:

```math
\\mathcal{C} = \\{ A_1 A_2 \\dots A_n ~|~ A_i \\in \\mathcal{A}_i \\}

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

MatrixZonotopeProduct(ms::MatrixZonotope{N,S}...) where {N,S} = MatrixZonotopeProduct(collect(ms))

"""
    *(A::MatrixZonotope{N,S}, B::MatrixZonotope{N,S}) where {N,S}

Appends matrix zonotope B to an existing matrix zonotope product P.
"""
function Base.:*(A::MatrixZonotope{N,S}, B::MatrixZonotope{N,S}) where {N,S}
    return MatrixZonotopeProduct([A, B])
end

"""
    *(P::MatrixZonotopeProduct{N,S}, B::MatrixZonotope{N,S}) where {N,S}

Appends matrix zonotope B to an existing matrix zonotope product P.
"""
function Base.:*(P::MatrixZonotopeProduct{N,S}, B::MatrixZonotope{N,S}) where {N,S}
    return MatrixZonotopeProduct(vcat(P.factors, B))
end

"""
    *(A::MatrixZonotope{N,S}, P::MatrixZonotopeProduct{N,S}) where {N,S}

Prepends matrix zonotope A to an existing matrix zonotope product P.
"""
function Base.:*(A::MatrixZonotope{N,S}, P::MatrixZonotopeProduct{N,S}) where {N,S}
    return MatrixZonotopeProduct(vcat(A, P.factors))
end

"""
    *(P1::MatrixZonotopeProduct{N,S}, P2::MatrixZonotopeProduct{N,S}) where {N,S}

Concatenates two matrix zonotope products.
"""
function Base.:*(P1::MatrixZonotopeProduct{N,S}, P2::MatrixZonotopeProduct{N,S}) where {N,S}
    return MatrixZonotopeProduct(vcat(P1.factors, P2.factors))
end

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

A matrix zonotope with redundant constant factors removed.
"""
function remove_redundant_factors(MZP::MatrixZonotopeProduct)
    gens = Vector{typeof(factors(MZP)[1])}()
    
    # right multiply the first n-1 factors
    @inbounds for (i, MZ) in enumerate(factors(MZP)[1:(end - 1)])
        if isempty(generators(MZ))
            G = linear_map(center(MZ), factors(MZP)[i + 1])
            gens[i] = G
        else
            push!(gens, MZ)
        end
    end

    # left multiply the last factor
    last_factor = factors(MZP)[end]
    if isempty(generators(last_factor))
        gens[end] = linear_map(MZ, center(last_factor))
    else
        push!(gens, last_factor)
    end

    return MatrixZonotopeProduct(gens)
end
