export SimpleSparsePolynomialZonotope, PolynomialZonotope, expmat, nparams,
       linear_map, quadratic_map, remove_redundant_generators

"""
    SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}} <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope that is *simple* in the sense that there is no distinction between independent and dependent generators.

A simple sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is represented by the set
```math
\\mathcal{PZ} = \\left\\{x \\in \\mathbb{R}^n : x = c + \\sum_{i=1}^h \\left(\\prod_{k=1}^p \\alpha_k^{E_{k, i}} \\right) g_i,~~ \\alpha_k \\in [-1, 1]~~ \\forall i = 1,\\ldots,p \\right\\},
```
where ``c ∈ \\mathbb{R}^n`` is the offset vector (or center), ``G ∈ \\mathbb{R}^{n \\times h}`` is the generator matrix with columns ``g_i``
(each ``g_i`` is called a *generator*), and where ``E ∈ \\mathbb{N}^{p×h}_{≥0}`` is the exponent matrix with matrix elements ``E_{k, i}``.

### Fields

- `c` -- offset vector
- `G` -- generator matrix
- `E` -- exponent matrix

### Notes

Sparse polynomial zonotopes were introduced in [KA21]. The *simple* variation was defined in [K22].

- [KA21] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*. Transactions on Automatic Control, 2021.
- [K21] N. Kochdumper. *Challenge Problem 5: Polynomial Zonotopes in Julia.* JuliaReach and JuliaIntervals Days 3, 2021.
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

Return the dimension of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The ambient dimension of `P`.
"""
dim(P::SSPZ) = size(P.c, 1)

"""
    ngens(P::SimpleSparsePolynomialZonotope)

Return the number of generators of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of generators of `P`.

### Notes

This is equivalent to the number of monomials in the polynomial representation of `P`.
"""
ngens(P::SSPZ) = size(P.G, 2)

"""
    nparams(P::SimpleSparsePolynomialZonotope)

Return the number of parameters in the polynomial representation of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of parameters in the polynomial representation of P.

### Notes

This corresponds to the number rows in the exponent matrix ``E`` (`p` in the set definition).

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

Return the order of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The order of `P`, defined as the quotient between the number of generators and the ambient dimension.
"""
order(P::SSPZ) = ngens(P) // dim(P)

"""
    center(P::SimpleSparsePolynomialZonotope)

Return the center of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The center of `P`.
"""
center(P::SSPZ) = P.c

"""
    genmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of generators of `P`.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of generators of `P`.
"""
genmat(P::SSPZ) = P.G

"""
    expmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of exponents of the sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``\alpha_k`` in the definition) and each column to a monomial.

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

Apply the linear map `M` to the simple sparse polynomial zonotope `P`.

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

Return an overapproximation of the quadratic map of the given polynomial zonotope.

### Input

- `Q` -- vector of square matrices
- `S` -- simple sparse polynomial zonotope

### Output

The quadratic map of the given zonotope represented as a polynomial zonotope.

### Algorithm

This method implements Proposition 12 in [1].
See also Proposition 3.1.30 in [2].

[1] N. Kochdumper, M. Althoff. *Sparse polynomial zonotopes: A novel set
representation for reachability analysis*. 2021
[2] N. Kochdumper. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems*. 2021.
"""
function quadratic_map(Q::Vector{MT}, S::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
    m = length(Q)
    c = center(S)
    h = ngens(S)
    G = genmat(S)
    E = expmat(S)

    cnew = similar(c, m)
    Gnew = similar(G, m, h^2 + h)
    QiG = similar(Q)
    @inbounds for (i, Qi) in enumerate(Q)
        cnew[i] = dot(c, Qi, c)
        Gnew[i, 1:h] = c' * (Qi + Qi') * G
        QiG[i] = Qi * G
    end

    Enew = repeat(E, 1, h + 1)
    @inbounds for i in 1:h
        idxstart = h * i + 1
        idxend = (i + 1) * h
        Enew[:, idxstart:idxend] .+= E[:, i]
        for j in eachindex(QiG)
            Gnew[j, idxstart:idxend] = G[:, i]' * QiG[j]
        end
    end
    return SimpleSparsePolynomialZonotope(cnew, Gnew, Enew)
end

"""
    remove_redundant_generators(S::SimpleSparsePolynomialZonotope)

Remove redundant generators from `S`.

### Input

- `S` -- simple sparse polynomial zonotope

### Output

A new simple simple sparse polynomial zonotope such that redundant generators have been reduced.

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the generator matrix and `E` the exponent matrix of `S`. The following simplifications are performed:

- Zero columns in `G` and the corresponding columns in `E` are eliminated.
- For zero columns in `E`, the corresponding column in `G` is summed to the center.
- Repeated columns in `E` are grouped together by summing the corresponding columns in `G`.
"""
function remove_redundant_generators(S::SimpleSparsePolynomialZonotope)

    c = center(S)
    G = genmat(S)
    E = expmat(S)

    Gnew = Matrix{eltype(G)}(undef, size(G, 1), 0)
    Enew = Matrix{eltype(E)}(undef, size(E, 1), 0)
    cnew = copy(c)

    visited_exps = Dict{Vector{Int}, Int}()
    @inbounds for (gi, ei) in zip(eachcol(G), eachcol(E))
        iszero(gi) && continue
        if iszero(ei)
            cnew += gi
        elseif haskey(visited_exps, ei) # repeated exponent
            idx = visited_exps[ei]
            Gnew[:, idx] += gi
        else
            Gnew = hcat(Gnew, gi)
            Enew = hcat(Enew, ei)
            visited_exps[ei] = size(Enew, 2)
        end
    end

    return SimpleSparsePolynomialZonotope(cnew, Gnew, Enew)
end
