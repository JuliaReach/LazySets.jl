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
"""
    overapproximate(P::SSPZ, ::Type{Zonotope}; nsdiv=1)

Returns a zonotope containing ``P``.

### Input

- `P`     -- simple sparse polynomial zonotope
- `nsdiv` -- (optional, default: `1`) size of uniform partitioning grid

### Output

A zonotope containing `P`.

"""
function overapproximate(P::SSPZ, ::Type{Zonotope}; nsdiv=1)
    (nsdiv != 1) && return overapproximate(P, UnionSetArray{Zonotope}; nsdiv=nsdiv)

    G = genmat(P)
    E = expmat(P)
    cnew = copy(center(P))
    Gnew = copy(G)
    @inbounds for (j, g) in enumerate(eachcol(G))
        if all(iseven, E[:, j])
            cnew .+= 0.5 * g
            Gnew[:, j] ./= 2
        end
    end
    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(P::SSPZ, ::Type{Zonotope}, dom::IntervalBox)

Compute the zonotope overapproximation of the given sparse polynomial zonotope over
the parameter domain `dom`, which should be a subset of `[-1, 1]^q`, where `q = nparams(P)`.
"""
function overapproximate(P::SSPZ, ::Type{Zonotope}, dom::IntervalBox)
    @assert dom ⊆ IntervalBox(-1..1, nparams(P)) "dom should be a subset of [-1, 1]^q"
    G = genmat(P)
    E = expmat(P)
    cnew = copy(center(P))
    Gnew = similar(G)
    @inbounds for (j, g) in enumerate(eachcol(G))
        α = IA.Interval(1, 1)
        for (i, vi) in enumerate(dom)
            α *= _fast_interval_pow(vi, E[i, j])
        end
        # α = mapreduce(x -> _fast_interval_pow(x[1],  x[2]), *, zip(dom, E[:, i])) # monomial value over the domain
        m, r = midpoint_radius(α)
        cnew .+= m * g
        Gnew[:, j] .= abs.(g) * r
    end
    return Zonotope(cnew, Gnew)
end

function overapproximate(P::SSPZ, ::Type{UnionSetArray{Zonotope}}; nsdiv=100)
    q = size(P.E, 1)
    dom = IntervalBox(IA.Interval(-1, 1), q)
    cells = mince(dom, nsdiv)
    return UnionSetArray([overapproximate(P, Zonotope, c) for c in cells])
end
