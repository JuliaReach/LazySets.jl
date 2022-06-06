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

"""
    nparams(P::SimpleSparsePolynomialZonotope)

Returns the number of variables used in the polynomial representation, which corresponds to
the number of rows in the exponents matrix.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of variables in the polynomial representation of P.

### Examples

```jldoctest
julia> S = SimpleSparsePolynomialZonotope([2.0, 0], [1 2;2 2.], [1 4;1 2])
SimpleSparsePolynomialZonotope{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Int64}}([2.0, 0.0], [1.0 2.0; 2.0 2.0], [1 4; 1 2])

julia> nparams(S)
2
```
"""
nparams(P::SSPZ) = size(P.E, 1)

order(P::SSPZ) = ngens(P) // dim(P)

center(P::SSPZ) = P.c

genmat(P::SSPZ) = P.G

"""
    expmat(P::SimpleSparsePolynomialZonotope)

Returns the matrix of exponents of the sparse polynomial zonotope. In the matrix, each
row corresponds to a variable and each column to a monomial.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Examples

```jldoctest
julia> S = SimpleSparsePolynomialZonotope([2.0, 0], [1 2;2 2.], [1 4;1 2])
SimpleSparsePolynomialZonotope{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Int64}}([2.0, 0.0], [1.0 2.0; 2.0 2.0], [1 4; 1 2])

julia> expmat(S)
2×2 Matrix{Int64}:
 1  4
 1  2
```
"""
expmat(P::SSPZ) = P.E
