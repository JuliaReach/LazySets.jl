export SimpleSparsePolynomialZonotope, PolynomialZonotope, expmat, nparams,
       linear_map, quadratic_map, remove_redundant_generators

"""
    SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N},
                                   MN<:AbstractMatrix{N},
                                   ME<:AbstractMatrix{<:Integer}}
        <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope that is *simple* in the sense
that there is no distinction between independent and dependent generators.

A simple sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n`` is
represented by the set
```math
\\mathcal{PZ} = \\left\\{x \\in \\mathbb{R}^n : x = c + \\sum_{i=1}^h \\left(\\prod_{k=1}^p \\alpha_k^{E_{k, i}} \\right) g_i,~~ \\alpha_k \\in [-1, 1]~~ \\forall i = 1,\\ldots,p \\right\\},
```
where ``c ∈ \\mathbb{R}^n`` is the offset vector (or center),
``G ∈ \\mathbb{R}^{n \\times h}`` is the generator matrix with columns ``g_i``
(each ``g_i`` is called a *generator*), and where ``E ∈ \\mathbb{N}^{p×h}_{≥0}``
is the exponent matrix with matrix elements ``E_{k, i}``.

### Fields

- `c` -- offset vector
- `G` -- generator matrix
- `E` -- exponent matrix

### Notes

Sparse polynomial zonotopes were introduced in [1]. The *simple* variation
was defined in [2].

- [1] N. Kochdumper and M. Althoff. *Sparse Polynomial Zonotopes: A Novel Set
Representation for Reachability Analysis*. Transactions on Automatic Control,
2021.
- [2] N. Kochdumper. *Challenge Problem 5: Polynomial Zonotopes in Julia.*
JuliaReach and JuliaIntervals Days 3, 2021.
"""
struct SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N},
        MN<:AbstractMatrix{N},
        ME<:AbstractMatrix{<:Integer}} <: AbstractPolynomialZonotope{N}
    c::VN
    G::MN
    E::ME

    function SimpleSparsePolynomialZonotope(c::VN, G::MN, E::ME) where {N,
            VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G " *
            "should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E " *
            "should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should contain " *
            "non-negative integers"))

        return new{N, VN, MN, ME}(c, G, E)
    end
end

"""
    PolynomialZonotope = SimpleSparsePolynomialZonotope

Alias for `SimpleSparsePolynomialZonotope`.

### Notes

Another shorthand is `SSPZ`.
"""
const PolynomialZonotope = SimpleSparsePolynomialZonotope

const SSPZ = SimpleSparsePolynomialZonotope

function isoperationtype(P::Type{<:SimpleSparsePolynomialZonotope})
    return false
end

function isconvextype(P::Type{<:SimpleSparsePolynomialZonotope})
    return false
end

function isboundedtype(P::Type{<:SimpleSparsePolynomialZonotope})
    return true
end

"""
    dim(P::SimpleSparsePolynomialZonotope)

Return the dimension of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The ambient dimension of `P`.
"""
dim(P::SSPZ) = length(P.c)

"""
    ngens(P::SimpleSparsePolynomialZonotope)

Return the number of generators of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of generators of `P`.

### Notes

This number corresponds to the number of monomials in the polynomial
representation of `P`.
"""
ngens(P::SSPZ) = size(P.G, 2)

"""
    nparams(P::SimpleSparsePolynomialZonotope)

Return the number of parameters in the polynomial representation of a simple
sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of parameters in the polynomial representation of P.

### Notes

This number corresponds to the number of rows in the exponent matrix ``E`` (`p`
in the mathematical set definition).

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

Return the order of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The order of `P`, defined as the quotient between the number of generators and
the ambient dimension.
"""
order(P::SSPZ) = ngens(P) // dim(P)

"""
    center(P::SimpleSparsePolynomialZonotope)

Return the center of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The center of `P`.
"""
center(P::SSPZ) = P.c

"""
    genmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of generators of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of generators of `P`.
"""
genmat(P::SSPZ) = P.G

"""
    expmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of exponents of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``\alpha_k`` in the
mathematical set definition) and each column corresponds to a monomial.

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
    linear_map(M::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling},
               P::SimpleSparsePolynomialZonotope)

Apply the linear map `M` to a simple sparse polynomial zonotope.

### Input

- `M` -- matrix
- `P` -- simple sparse polynomial zonotope

### Output

The set resulting from applying the linear map `M` to `P`.
"""
function linear_map(M::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling},
                    P::SSPZ)
    return SimpleSparsePolynomialZonotope(M * center(P), M * genmat(P), expmat(P))
end

"""
    quadratic_map(Q::Vector{MT}, S::SimpleSparsePolynomialZonotope)
        where {N, MT<:AbstractMatrix{N}}

Return the quadratic map of a simple sparse polynomial zonotope.

### Input

- `Q` -- vector of square matrices
- `S` -- simple sparse polynomial zonotope

### Output

The quadratic map of `P` represented as a simple sparse polynomial zonotope.

### Algorithm

This method implements Proposition 12 in [1].
See also Proposition 3.1.30 in [2].

[1] N. Kochdumper, M. Althoff. *Sparse polynomial zonotopes: A novel set
representation for reachability analysis*. 2021
[2] N. Kochdumper. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems*. 2021.
"""
function quadratic_map(Q::Vector{MT},
                       S::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
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
    Z = SimpleSparsePolynomialZonotope(cnew, Gnew, Enew)
    return remove_redundant_generators(Z)
end


"""
    quadratic_map(Q::Vector{MT}, S1::SimpleSparsePolynomialZonotope,
                  S2::SimpleSparsePolynomialZonotope)
        where {N, MT<:AbstractMatrix{N}}

Return the quadratic map of two simple sparse polynomial zonotopes.
The quadratic map is the set
``\\{x | xᵢ = s₁ᵀQᵢs₂, s₁ ∈ S₁, s₂ ∈ S₂, Qᵢ ∈ Q\\}``.

### Input

- `Q`  -- vector of square matrices
- `S1` -- simple sparse polynomial zonotope
- `S2` -- simple sparse polynomial zonotope

### Output

The quadratic map of the given simple sparse polynomial zonotopes represented as
a simple sparse polynomial zonotope.

### Algorithm

This method implements Proposition 3.1.30 in [1].

[1] N. Kochdumper. *Extensions of polynomial zonotopes and their application to
verification of cyber-physical systems*. 2021.
"""
function quadratic_map(Q::Vector{MT}, S1::SimpleSparsePolynomialZonotope,
                       S2::SimpleSparsePolynomialZonotope
                      ) where {N, MT<:AbstractMatrix{N}}
    @assert nparams(S1) == nparams(S2)

    c1 = center(S1)
    c2 = center(S2)
    G1 = genmat(S1)
    G2 = genmat(S2)
    E1 = expmat(S1)
    E2 = expmat(S2)

    c = [dot(c1, Qi, c2) for Qi in Q]

    Ghat1 = reduce(vcat, c2' * Qi' * G1 for Qi in Q)
    Ghat2 = reduce(vcat, c1' * Qi * G2 for Qi in Q)

    Gbar = reduce(hcat, reduce(vcat, gj' * Qi * G2 for Qi in Q) for gj in eachcol(G1))
    Ebar = reduce(hcat, E2 .+ e1j for e1j in eachcol(E1))

    G = hcat(Ghat1, Ghat2, Gbar)
    E = hcat(E1, E2, Ebar)

    return remove_redundant_generators(SimpleSparsePolynomialZonotope(c, G, E))
end

"""
    remove_redundant_generators(S::SimpleSparsePolynomialZonotope)

Remove redundant generators from a simple sparse polynomial zonotope.

### Input

- `S` -- simple sparse polynomial zonotope

### Output

A new simple sparse polynomial zonotope such that redundant generators have been
removed.

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the generator matrix and `E` the exponent matrix of `S`. The
following simplifications are performed:

- Zero columns in `G` and the corresponding columns in `E` are removed.
- For zero columns in `E`, the corresponding column in `G` is summed to the
  center.
- Repeated columns in `E` are grouped together by summing the corresponding
  columns in `G`.
"""
function remove_redundant_generators(S::SimpleSparsePolynomialZonotope)

    c, G, E = _remove_redundant_generators_polyzono(center(S), genmat(S), expmat(S))

    return SimpleSparsePolynomialZonotope(c, G, E)
end

function _remove_redundant_generators_polyzono(c, G, E)
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

    return cnew, Gnew, Enew
end

"""
    rand(::Type{SimpleSparsePolynomialZonotope};
         [N]::Type{<:Real}=Float64, [dim]::Int=2, [nparams]::Int=2,
         [maxdeg]::Int=3, [num_generators]::Int=-1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random simple sparse polynomial zonotope.

### Input

- `Zonotope`       -- type for dispatch
- `N`              -- (optional, default: `Float64`) numeric type
- `dim`            -- (optional, default: 2) dimension
- `nparams`        -- (optional, default: 2) number of parameters
- `maxdeg`         -- (optinal, default: 3) maximum degree for each parameter
- `rng`            -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`           -- (optional, default: `nothing`) seed for reseeding
- `num_generators` -- (optional, default: `-1`) number of generators of the
                      zonotope (see comment below)

### Output

A random simple sparse polynomial zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the argument `num_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator). Note that the
final number of generators may be lower if redundant monomials are generated.
"""
function rand(::Type{SimpleSparsePolynomialZonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              nparams::Int=2,
              maxdeg::Int=3,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_generators::Int=-1)
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    if num_generators < 0
        num_generators = (dim == 1) ? 1 : rand(rng, dim:2*dim)
    end
    generators = randn(rng, N, dim, num_generators)
    expmat = rand(rng, 0:maxdeg, nparams, num_generators)
    SSPZ = SimpleSparsePolynomialZonotope(center, generators, expmat)
    return remove_redundant_generators(SSPZ)
end

"""
    convex_hull(P::SimpleSparsePolynomialZonotope)

Compute the convex hull of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The tightest convex simple sparse polynomial zonotope containing `P`.
"""
function convex_hull(P::SimpleSparsePolynomialZonotope)
    return linear_combination(P, P)
end
