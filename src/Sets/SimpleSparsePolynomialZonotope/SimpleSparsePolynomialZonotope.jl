"""
    SimpleSparsePolynomialZonotope{N, VN<:AbstractVector{N},
                                   MN<:AbstractMatrix{N},
                                   ME<:AbstractMatrix{Int}}
        <: AbstractSparsePolynomialZonotope{N}

Type that represents a sparse polynomial zonotope that is *simple* in the sense
that there is no distinction between independent and dependent generators.

A simple sparse polynomial zonotope ``\\mathcal{PZ} ⊂ ℝ^n`` is
represented by the set
```math
\\mathcal{PZ} = \\left\\{x ∈ ℝ^n : x = c + ∑_{i=1}^h \\left(∏_{k=1}^p α_k^{E_{k, i}} \\right) g_i,~~ α_k ∈ [-1, 1]~~ ∀ i = 1,…,p \\right\\},
```
where ``c ∈ ℝ^n`` is the offset vector (or center),
``G ∈ ℝ^{n × h}`` is the generator matrix with columns ``g_i``
(each ``g_i`` is called a *generator*), and where ``E ∈ \\mathbb{N}^{p×h}_{≥0}``
is the exponent matrix with matrix elements ``E_{k, i}``.

### Fields

- `c` -- offset vector
- `G` -- generator matrix
- `E` -- exponent matrix

### Notes

Sparse polynomial zonotopes were introduced in [KochdumperA21](@citet).
The *simple* variant was defined in [Kochdumper21b].
"""
struct SimpleSparsePolynomialZonotope{N,VN<:AbstractVector{N},
                                      MN<:AbstractMatrix{N},
                                      ME<:AbstractMatrix{Int}} <:
       AbstractSparsePolynomialZonotope{N}
    c::VN
    G::MN
    E::ME

    function SimpleSparsePolynomialZonotope(c::VN, G::MN,
                                            E::ME) where {N,
                                                          VN<:AbstractVector{N},
                                                          MN<:AbstractMatrix{N},
                                                          ME<:AbstractMatrix{Int}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G " *
                                                                "should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E " *
                                                                 "should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should contain " *
                                                  "non-negative integers"))

        return new{N,VN,MN,ME}(c, G, E)
    end
end
