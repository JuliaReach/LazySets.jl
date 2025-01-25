"""
    Star{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}} <: AbstractPolyhedron{N}

Generalized star set with a polyhedral predicate, i.e.

```math
X = \\{x ∈ ℝ^n : x = x₀ + ∑_{i=1}^m α_i v_i,~~\\textrm{s.t. } P(α) = ⊤ \\},
```
where ``x₀ ∈ ℝ^n`` is the center, the ``m`` vectors ``v₁, …, vₘ`` form
the basis of the star set, and the combination factors
``α = (α₁, …, αₘ) ∈ ℝ^m`` are the predicate's decision variables,
i.e., ``P : α ∈ ℝ^m → \\{⊤, ⊥\\}`` where the polyhedral predicate
satisfies ``P(α) = ⊤`` if and only if ``A·α ≤ b`` for some fixed
``A ∈ ℝ^{p × m}`` and ``b ∈ ℝ^p``.

### Fields

- `c` -- vector that represents the center
- `V` -- matrix where each column corresponds to a basis vector
- `P` -- polyhedral set that represents the predicate

### Notes

The predicate function is implemented as a conjunction of linear constraints,
i.e., a subtype of `AbstractPolyhedron`. By a slight abuse of notation, the
predicate is also used to denote the subset of ``ℝ^n`` such that
``P(α) = ⊤`` holds.

The ``m`` basis vectors (each one ``n``-dimensional) are stored as the columns
of an ``n × m`` matrix.

We remark that a `Star` is mathematically equivalent to the affine map of the
polyhedral set `P`, with the transformation matrix and translation vector being
`V` and `c`, respectively.

Star sets as defined here were introduced in [DuggiralaV16](@citet); see also
[BakD17](@citet) for a preliminary definition.

### Examples

This example is drawn from Example 1 in [2]. Consider the two-dimensional plane
``ℝ^2``. Let

```jldoctest star_constructor
julia> V = [[1.0, 0.0], [0.0, 1.0]];
```
be the basis vectors and take

```jldoctest star_constructor
julia> c = [3.0, 3.0];
```
as the center of the star set. Let the predicate be the infinity-norm ball of
radius 1,

```jldoctest star_constructor
julia> P = BallInf(zeros(2), 1.0);
```
We construct the star set ``X = ⟨c, V, P⟩`` as follows:

```jldoctest star_constructor
julia> S = Star(c, V, P)
Star{Float64, Vector{Float64}, Matrix{Float64}, BallInf{Float64, Vector{Float64}}}([3.0, 3.0], [1.0 0.0; 0.0 1.0], BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0))
```

We can use getter functions for each component field:

```jldoctest star_constructor
julia> center(S)
2-element Vector{Float64}:
 3.0
 3.0

julia> basis(S)
2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> predicate(S)
BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0)
```
In this case, the generalized star ``S`` above is equivalent to the rectangle
``T`` below.

```math
    T = \\{(x, y) ∈ ℝ^2 : (2 ≤ x ≤ 4) ∧ (2 ≤ y ≤ 4) \\}
```
"""
struct Star{N,VN<:AbstractVector{N},MN<:AbstractMatrix{N},PT<:AbstractPolyhedron{N}} <:
       AbstractPolyhedron{N}
    c::VN # center
    V::MN # basis
    P::PT # predicate

    # default constructor with size checks
    function Star(c::VN, V::MN,
                  P::PT) where {N,VN<:AbstractVector{N},MN<:AbstractMatrix{N},
                                PT<:AbstractPolyhedron{N}}
        @assert length(c) == size(V, 1) "a center of length $(length(c)) is " *
                                        "incompatible with basis vectors of length $(size(V, 1))"

        @assert dim(P) == size(V, 2) "the number of basis vectors " *
                                     "$(size(V, 2)) is incompatible with the predicate's dimension " *
                                     "$(dim(P))"

        return new{N,VN,MN,PT}(c, V, P)
    end
end

# constructor from list of generators
function Star(c::VN, vlist::AbstractVector{VN},
              P::PT) where {N,VN<:AbstractVector{N},PT<:AbstractPolyhedron{N}}
    V = to_matrix(vlist, length(c))
    return Star(c, V, P)
end
