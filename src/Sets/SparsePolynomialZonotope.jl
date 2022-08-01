export SparsePolynomialZonotope, expmat, nparams, ndependentgens, nindependentgens,
       dependent_genmat, independent_genmat, indexvector,
       linear_map, quadratic_map, remove_redundant_generators

"""
    SparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, MNI<:AbstractVecOrMat{N}, ME<:AbstractMatrix{<:Integer}, VI<:AbstractVector{<:Integer}} <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope.

A sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is represented by the set
```math
\\mathcal{PZ} = \\left\\{x \\in \\mathbb{R}^n : x = c + ∑ᵢ₌₁ʰ\\left(∏ₖ₌₁ᵖ α_k^{E_{k, i}} \\right)Gᵢ+∑ⱼ₌₁^qβⱼGIⱼ,~~ α_k ∈ [-1, 1]~~ ∀ i = 1,…,p, j=1,…,q \\right\\},
```
where ``c ∈ \\mathbb{R}^n`` is the offset vector (or center), ``G ∈ \\mathbb{R}^{n \\times h}`` is the dependent generator matrix with columns ``Gᵢ``,
``GI ∈ \\mathbb{R}^{n×q}`` is the independent generator matrix and ``E ∈ \\mathbb{N}^{p×h}_{≥0}`` is the exponent matrix with matrix elements ``E_{k, i}``.

### Fields

- `c`   -- offset vector
- `G`   -- dependent generator matrix
- `GI`  -- independent generator matrix
- `E`   -- exponent matrix
- `idx` -- identifier vector, vector of positive integers identifing the dependent parameters of `PZ`.

### Notes

Sparse polynomial zonotopes were introduced in [KA21].

- [KA21] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*. Transactions on Automatic Control, 2021.
"""
struct SparsePolynomialZonotope{N,
                                VN<:AbstractVector{N},
                                MN<:AbstractMatrix{N},
                                MNI<:AbstractVecOrMat{N},
                                ME<:AbstractMatrix{<:Integer},
                                VI<:AbstractVector{<:Integer}} <: AbstractPolynomialZonotope{N}
    c::VN
    G::MN
    GI::MNI
    E::ME
    idx::VI

    function SparsePolynomialZonotope(c::VN, G::MN, GI::MNI, E::ME, idx::VI) where {N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, MNI<:AbstractVecOrMat{N}, ME<:AbstractMatrix{<:Integer}, VI<:AbstractVector{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G should have the same number of rows"))
        @assert length(c) == size(GI, 1) throw(DimensionMismatch("c and GI should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should have non-negative integers"))
        @assert all(>(0), idx) throw(ArgumentError("identifiers in index vector must be positive integers"))
        return new{N, VN, MN, MNI, ME, VI}(c, G, GI, E, idx)
    end
end

function SparsePolynomialZonotope(c::AbstractVector, G::AbstractMatrix, GI::AbstractVecOrMat, E::AbstractMatrix{<:Integer})
    n = size(E, 1)
    return SparsePolynomialZonotope(c, G, GI, E, 1:n)
end


const SPZ = SparsePolynomialZonotope

"""
    dim(P::SparsePolynomialZonotope)

Return the dimension of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The ambient dimension of `P`.
"""
dim(P::SPZ) = length(P.c)

"""
    ndependentgens(P::SparsePolynomialZonotope)

Return the number of dependent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
ndependentgens(P::SPZ) = size(P.G, 2)

"""
    nindependentgens(P::SparsePolynomialZonotope)

Return the number of independent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
nindependentgens(P::SPZ) = size(P.GI, 2)

"""
    nparams(P::SparsePolynomialZonotope)

Return the number of dependent parameters in the polynomial representation of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The number of dependent parameters in the polynomial representation of P.

### Notes

This corresponds to the number rows in the exponent matrix ``E``.
"""
nparams(P::SPZ) = size(P.E, 1)

"""
    order(P::SparsePolynomialZonotope)

Return the order of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The order of `P`, defined as the quotient between the number of generators and the ambient dimension.
"""
order(P::SPZ) = (ndependentgens(P) + nindependentgens(P)) // dim(P)

"""
    center(P::SparsePolynomialZonotope)

Return the center of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Output

The center of `P`.
"""
center(P::SPZ) = P.c

"""
    dependent_genmat(P::SparsePolynomialZonotope)

Return the matrix of dependent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
dependent_genmat(P::SPZ) = P.G

"""
    independent_genmat(P::SparsePolynomialZonotope)

Return the matrix of independent generators of `P`.

### Input

- `P` -- sparse polynomial zonotope
"""
independent_genmat(P::SPZ) = P.GI

"""
    expmat(P::SparsePolynomialZonotope)

Return the matrix of exponents of the sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``\alpha_k`` in the definition) and each column to a monomial.
"""
expmat(P::SPZ) = P.E

"""
    indexvector(P::SparsePolynomialZonotope)

Return the index vector of `P`.

### Input

- `P` -- sparse polynomial zonotope

### Notes

The index vector is a vector of positive integers identifing the dependent parameters of `P`.
"""
indexvector(P::SPZ) = P.idx

"""
    linear_map(M::AbstractMatrix, P::SparsePolynomialZonotope)

Apply the linear map `M` to the sparse polynomial zonotope `P`.

### Input

- `M` -- square matrix with `size(M) == dim(P)`
- `P` -- sparse polynomial zonotope

### Output

The set resulting from applying the linear map `M` to `P`.
"""
function linear_map(M::AbstractMatrix, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * dependent_genmat(P),
                                    M * independent_genmat(P),
                                    expmat(P),
                                    indexvector(P)
                                    )
end

# function quadratic_map(Q::Vector{MT}, S::SparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
#     m = length(Q)
#     c = center(S)
#     h = ngens(S)
#     G = genmat(S)
#     E = expmat(S)

#     cnew = similar(c, m)
#     Gnew = similar(G, m, h^2 + h)
#     QiG = similar(Q)
#     @inbounds for (i, Qi) in enumerate(Q)
#         cnew[i] = dot(c, Qi, c)
#         Gnew[i, 1:h] = c' * (Qi + Qi') * G
#         QiG[i] = Qi * G
#     end

#     Enew = repeat(E, 1, h + 1)
#     @inbounds for i in 1:h
#         idxstart = h * i + 1
#         idxend = (i + 1) * h
#         Enew[:, idxstart:idxend] .+= E[:, i]
#         for j in eachindex(QiG)
#             Gnew[j, idxstart:idxend] = G[:, i]' * QiG[j]
#         end
#     end
#     return SparsePolynomialZonotope(cnew, Gnew, Enew)
# end

# function remove_redundant_generators(S::SparsePolynomialZonotope)

#     c = center(S)
#     G = genmat(S)
#     E = expmat(S)

#     Gnew = Matrix{eltype(G)}(undef, size(G, 1), 0)
#     Enew = Matrix{eltype(E)}(undef, size(E, 1), 0)
#     cnew = copy(c)

#     visited_exps = Dict{Vector{Int}, Int}()
#     @inbounds for (gi, ei) in zip(eachcol(G), eachcol(E))
#         iszero(gi) && continue
#         if iszero(ei)
#             cnew += gi
#         elseif haskey(visited_exps, ei) # repeated exponent
#             idx = visited_exps[ei]
#             Gnew[:, idx] += gi
#         else
#             Gnew = hcat(Gnew, gi)
#             Enew = hcat(Enew, ei)
#             visited_exps[ei] = size(Enew, 2)
#         end
#     end

#     return SparsePolynomialZonotope(cnew, Gnew, Enew)
# end
