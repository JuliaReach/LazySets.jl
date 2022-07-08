export SimpleSparsePolynomialZonotope, PolynomialZonotope, expmat, nparams,
       linear_map, quadratic_map

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

"""
    PolynomialZonotope = SimpleSparsePolynomialZonotope

Alias for `SimpleSparsePolynomialZonotope`.
"""
const PolynomialZonotope = SimpleSparsePolynomialZonotope


const SSPZ = SimpleSparsePolynomialZonotope

"""
    dim(P::SimpleSparsePolynomialZonotope)

Computes the dimension of the simple sparse polynomial zonotope `P`.
"""
dim(P::SSPZ) = size(P.c, 1)

"""
    ngens(P::SimpleSparsePolynomialZonotope)

Computes the number of generators of the simple sparse polynomial zonotope `P`. This is
equivalent to the number of monomials in the polynomial representation of `P`.
"""
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

"""
    order(P::SimpleSparsePolynomialZonotope)

Computes the order of the simple sparse polynomial zonotope `P`.
"""
order(P::SSPZ) = ngens(P) // dim(P)

"""
    center(P::SimpleSparsePolynomialZonotope)

Returns the center of the simple sparse polynomial zonotope `P`.
"""
center(P::SSPZ) = P.c

"""
    genmat(P::SimpleSparsePolynomialZonotope)

Returns the matrix of generators of the simple sparse polynomial zonotope `P`.
"""
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

"""
    linear_map(M::AbstractMatrix, P::SimpleSparsePolynomialZonotope)

applies the linear mapping `M` to the simple sparse polynomial zonotope `P`.

### Input

- `M` -- square matrix with size(M) == dim(P)
- `P` -- simple sparse polynomial zonotope

### Output

The set resulting from applying the linear map `M` to `P`.
"""
function linear_map(M::AbstractMatrix, P::SSPZ)
    return SimpleSparsePolynomialZonotope(M * center(P), M * genmat(P), expmat(P))
end

"""
    quadratic_map(Q::Vector{MT}, S::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}

Return an overapproximation of the quadratic map of the given zonotope.

### Input

- `S` -- simple sparse polynomial zonotope
- `Q` -- array of square matrices

### Output

An overapproximation of the quadratic map of the given zonotope.
"""
function quadratic_map(Q::Vector{MT}, S::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
    m = length(Q)
    c = center(S)
    h = ngens(S)
    G = genmat(S)
    E = expmat(S)

    cnew = similar(c, m)
    Gnew = similar(G, m, h^2 + h)
    tmp = similar(Q)
    @inbounds for (i, q) in enumerate(Q)
        cnew[i] = dot(c, q, c)
        Gnew[i, 1:h] = c' * (q + q')*G
        tmp[i] = q * G
    end

    Enew = repeat(E, 1, h + 1)
    @inbounds for i in 1:h
        idxstart = h * i + 1
        idxend = (i + 1) * h
        Enew[:, idxstart:idxend] .+= E[:, i]
        for j in eachindex(tmp)
            Gnew[j, idxstart:idxend] = G[:, i]' * tmp[j]
        end
    end
    return SimpleSparsePolynomialZonotope(cnew, Gnew, Enew)
end
