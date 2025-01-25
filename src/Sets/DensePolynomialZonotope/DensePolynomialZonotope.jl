"""
    DensePolynomialZonotope{N, VT, VMT, MT} <: AbstractPolynomialZonotope{N}

Type that represents a polynomial zonotope.

### Fields

- `c`  -- starting point
- `E`  -- matrix of multi-indexed generators such that *all* indices have the
          same value
- `F`  -- matrix of multi-indexed generators such that *not all* indices have
          the same value
- `G`  -- matrix of single-indexed generators

### Notes

Polynomial zonotopes were introduced in [Althoff13](@citet) and have been applied as non-convex
set representation to the reachability problem of nonlinear ODEs.

Mathematically, a polynomial zonotope is defined as the tuple ``(c, E, F, G)``,
where:

- ``c ∈ ℝ^n`` is the starting point (in some particular cases it
  corresponds to the center of a usual zonotope),

- ``E = [E^{[1]} ⋯ E^{[η]}]`` is an ``n × p × η`` matrix with column
  blocks

```math
E^{[i]} = [f^{([i], 1, 1, …, 1)} ⋯ f^{([i], p, p, …, p)}], \\qquad i = 1,…, η
```
called the matrix of *multi-indexed generators with equal indices*, where each
``f^{([i], k_1, k_2, …, k_i)}`` is an ``n``-vector,

- ``F = [F^{[2]} ⋯ F^{[η]}]`` is a matrix with column-blocks

```math
F^{[i]} = [f^{([i], 1, 1, …, 1, 2)} f^{([i], 1, 1, …, 1, 3)} ⋯ f^{([i], 1, 1, …, 1, p)} \\\\
f^{([i], 1, 1, …, 2, 2)} f^{([i], 1, 1, …, 2, 3)} ⋯ f^{([i], 1, 1, …, 2, p)} \\\\
f^{([i], 1, 1, …, 3, 3)} ⋯], \\qquad i = 1,…, η
```
called the matrix of *multi-indexed generators with unequal indices* (or, more
accurately, not-all-equal indices), where each ``f^{([i], k_1, k_2, …, k_i)}``
is an ``n``-vector,

- ``G = [G^{[1]} ⋯ G^{[q]}]`` is an ``n × q`` matrix with columns

```math
G^{[i]} = g^{(i)}, \\qquad i = 1,…, q
```
called the matrix of *single-indexed generators*, where each ``g^{(i)}`` is an
``n``-vector.

The polynomial zonotope ``(c, E, F, G)`` defines the set:

```math
\\left\\{ c + ∑_{j=1}^p β_j f^{([1], j)} + ∑_{j=1}^p ∑_{k=j}^p β_j β_k f^{([2], j, k)} + \\\\
+ … + ∑_{j=1}^p ∑_{k=j}^p ⋯ ∑_{m=ℓ}^p β_j β_k ⋯ β_m f^{([η], j, k, …, m)} + \\\\
+ ∑_{i=1}^q γ_i g^{(i)}, \\qquad β_i, γ_i ∈ [-1, 1] \\right\\},
```
where the number of factors in the final product, ``β_j β_k ⋯ β_m``, corresponds
to the polynomial order ``η``.
"""
struct DensePolynomialZonotope{N,VT,VMT,MT} <: AbstractPolynomialZonotope{N}
    c::VT
    E::VMT
    F::VMT
    G::MT

    # default constructor with dimension check
    function DensePolynomialZonotope(c::VT, E::VMT, F::VMT, G::MT) where {VT,VMT,MT}
        # check polynomial order
        η = length(E)
        @assert length(F) == (η > 0 ? η - 1 : 0) "incompatible lengths of E and F: $η and " *
                                                 "$(length(F))"
        # check dimension
        n = length(c)
        @assert size(G, 1) == n "incompatible dimension $n for G = $G"
        @assert all(size(A, 1) == n for A in E) "incompatible dimension $n for E = $E"

        N = eltype(c)
        return new{N,VT,VMT,MT}(c, E, F, G)
    end
end
