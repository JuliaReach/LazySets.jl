"""
    SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}} <: AbstractPolynomialZonotope{N}

Type that represents a sparse polynomial zonotope.

Given a constant offset ``c ∈ \\mathbb{R}^n``, a generator matrix ``G ∈ \\mathbb{R}^{n \times h}`` and an exponent matrix
``E ∈ \\mathbb{N}^{q×h}_{≥0}``, a sparse polynomial zonotope ``\\mathcal{PZ} ⊂ \\mathbb{R}^n``.

### Fields

- `c` -- constant offset vector
- `G` -- generators matrix
- `E` -- exponents matrix

### Notes

Sparse polynomial zonotopes were introduced in N. Kochdumper and M. Althoff.
*Sparse Polynomial Zonotopes: A Novel Set Representation for Reachability Analysis*. Transactions on Automatic Control, 2021.
"""
struct SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}}
    c::VN
    G::MN
    E::ME

    function SimpleSparsePolynomialZonotope(c::VN, G::MN, E::ME) where {N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}}
        length(c) == size(G, 1) || throw(DimensionMismatch("c and G should have the same number of rows"))
        size(G, 2) == size(E, 2) || throw(DimensionMismatch("G and E should have the same number of columns"))
        all(>=(0), E) || throw(ArgumentError("E should have non-negative integers"))

        return new{N, VN, MN, ME}(c, G, E)
    end
end

const SSPZ = SimpleSparsePolynomialZonotope

dim(p::SSPZ) = size(p.c, 1)
ngens(p::SSPZ) = size(p.G, 2)
nparams(p::SSPZ) = size(p.E, 1)
order(p::SSPZ) = ngens(p) // dim(p)

center(p::SSPZ) = p.c
genmat(p::SSPZ) = p.G
expmat(p::SSPZ) = p.E

"""
    overapproximate(p::SSPZ, ::Type{Zonotope}; nsdiv=1)

Returns a zonotope containing ``p``.

### Input

- `p`     -- simple sparse polynomial zonotope
- `nsdiv` -- (optional, default: `1`)

### Output

A zonotope containing `p`.

"""
function overapproximate(p::SSPZ, ::Type{Zonotope}; nsdiv=1)
    (nsdiv != 1) && return overapproximate(p, UnionSetArray{Zonotope}; nsdiv=nsdiv)

    G = genmat(p)
    E = expmat(p)
    cnew = copy(center(p))
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
    overapproximate(p::SSPZ, ::Type{Zonotope}, dom::IntervalBox)

Compute the zonotope overapproximation of the given sparse polynomial zonotope over
the parameter domain `dom`, which should be a subset of `[-1, 1]^q`, where `q = nparams(p)`.
"""
function overapproximate(p::SSPZ, ::Type{Zonotope}, dom::IntervalBox)
    G = genmat(p)
    E = expmat(p)
    cnew = copy(center(p))
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

function overapproximate(p::SSPZ, ::Type{UnionSetArray{Zonotope}}; nsdiv=100)
    q = size(p.E, 1)
    dom = IntervalBox(IA.Interval(-1, 1), q)
    cells = mince(dom, nsdiv)
    return UnionSetArray([overapproximate(p, Zonotope, c) for c in cells])
end
