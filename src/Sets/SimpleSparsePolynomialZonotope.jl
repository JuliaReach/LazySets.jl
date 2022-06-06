export SimpleSparsePolynomialZonotope, expmat, nparams

"""
    SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}} <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope.

A sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is represented by
a constant offset ``c ∈ \\mathbb{R}^n``, a generator matrix ``G ∈ \\mathbb{R}^{n \times h}`` and an exponent matrix
``E ∈ \\mathbb{N}^{q×h}_{≥0}``.

### Fields

- `c` -- constant offset vector
- `G` -- generator matrix
- `E` -- exponent matrix

### Notes

Sparse polynomial zonotopes were introduced in N. Kochdumper and M. Althoff.
*Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*. Transactions on Automatic Control, 2021.
"""
struct SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}} <: AbstractPolynomialZonotope{N}
    c::VN
    G::MN
    E::ME

    function SimpleSparsePolynomialZonotope(c::VN, G::MN, E::ME) where {N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should have non-negative integers"))

        return new{N, VN, MN, ME}(c, G, E)
    end
end

const SSPZ = SimpleSparsePolynomialZonotope

dim(P::SSPZ) = size(P.c, 1)
ngens(P::SSPZ) = size(P.G, 2)
nparams(P::SSPZ) = size(P.E, 1)
order(P::SSPZ) = ngens(P) // dim(P)

center(P::SSPZ) = P.c
genmat(P::SSPZ) = P.G
expmat(P::SSPZ) = P.E
