export Star,
       center,
       basis,
       predicate

"""
    Star{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}} <: AbstractStar{N}

Type that represents a generalized star set where the predicate is polyhedral, i.e.

```math
X = \\{x ∈ \\mathbb{R}^n : x = x₀ + \\sum_{i=1}^m α_i v_i,~~\\textrm{s.t. } P(α) = ⊤ \\},
```
where ``x₀ ∈ \\mathbb{R}^n`` is the center, the ``m`` vectors ``v₁, …, vₘ`` form
the basis of the star set, and the combination factors
``α = (α₁, …, αₘ) ∈ \\mathbb{R}^m`` are the predicates' decision variables,
i.e. ``P : α ∈ \\mathbb{R}^m → \\{⊤, ⊥\\}`` where the polyhedral predicate satisfies
``P(α) = ⊤`` if and only if ``Aα ≤ b`` for some fixed ``A ∈ \\mathbb{R}^{p × m}`` and
``b ∈ \\mathbb{R}^p``.

### Fields

- `c` -- vector that represents the center
- `V` -- matrix where each column corresponds to a basis vector
- `P` -- polyhedral set that represents the predicate

### Notes

The predicate function is implemented as a conjunction of linear constraints, i.e.
a subtype of `AbstractPolyhedron`. By a slight abuse of notation, the *predicate*
is also used to denote the subset of ``\\mathbb{R}^n`` such that ``P(α) = ⊤`` holds.

The ``m`` basis vectors (each one ``n``-dimensional) are stored as the columns
of an ``n × m`` matrix.

### Examples

This example is drawn from Example 1 in [2]. Consider the two-dimensional plane
``\\mathbb{R}^2``. Let

```jldoctest star_constructor
julia> V = [[1.0, 0.0], [0.0, 1.0]];
```
be the basis vectors and take

```jldoctest star_constructor
julia> c = [3.0, 3.0];
```
as the center of the star set. Let the predicate be the infinity-norm ball of radius 1,

```jldoctest star_constructor
julia> P = BallInf(zeros(2), 1.0);
```
Finally, the star set ``X = ⟨c, V, P⟩`` defines the set:

```jldoctest star_constructor
julia> S = Star(c, V, P)
Star{Float64,Array{Float64,1},Array{Float64,2},BallInf{Float64,Array{Float64,1}}}([3.0, 3.0], [1.0 0.0; 0.0 1.0], BallInf{Float64,Array{Float64,1}}([0.0, 0.0], 1.0))
```

We can use getter functions for each component field:

```jldoctest star_constructor
julia> center(S)
2-element Array{Float64,1}:
 3.0
 3.0

julia> basis(S)
2×2 Array{Float64,2}:
 1.0  0.0
 0.0  1.0

julia> predicate(S)
BallInf{Float64,Array{Float64,1}}([0.0, 0.0], 1.0)
```
In this case, we know calculating by hand that the generalized star ``S`` is
defined by the rectangular set

```math
    T = \\{(x, y) ∈ \\mathbb{R}^2 : (2 ≤ x ≤ 4) ∧ (2 ≤ y ≤ 4) \\}
```
It holds that ``T`` and ``S`` are equivalent.

### References

Star sets as defined here were introduced in [1]; see also [2] for a preliminary
definition. For applications in reachability analysis of neural networks, see
[3].

- [1] Duggirala, P. S., and Mahesh V. *Parsimonious, simulation based verification of linear systems.*
      International Conference on Computer Aided Verification. Springer, Cham, 2016.

- [2] Bak S, Duggirala PS. *Simulation-equivalent reachability of large linear systems with inputs.*
      In International Conference on Computer Aided Verification
      2017 Jul 24 (pp. 401-420). Springer, Cham.

- [3] Tran, H. D., Lopez, D. M., Musau, P., Yang, X., Nguyen, L. V., Xiang, W., & Johnson, T. T. (2019, October).
      *Star-based reachability analysis of deep neural networks.*
      In International Symposium on Formal Methods (pp. 670-686). Springer, Cham.
"""
struct Star{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}} <: AbstractStar{N}
    c::VN # center
    V::MN # basis
    P::PT # predicate

    # default constructor with size checks
    function Star(c::VN, V::MN, P::PT) where {N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}}
        @assert length(c) == size(V, 1) "the center of the basis vectors should be compatible, " *
                                        "but they are of length $(length(c)) and $(size(V, 1)) respectively"

        @assert dim(P) == size(V, 2) "the number of basis vectors should be compatible " *
                                     "with the predicates' dimension, but they are $(size(V, 2)) and $(dim(P)) respectively"

        return new{N, VN, MN, PT}(c, V, P)
    end
end

# constructor from center and list of generators
function Star(c::VN, Vlist::AbstractVector{VN}, P::PT) where {N, VN<:AbstractVector{N}, PT<:AbstractPolyhedron{N}}
    MT = _matrix_type(VN)
    V = MT(undef, length(c), length(Vlist))
    for (j, vj) in enumerate(Vlist)
        @inbounds V[:, j] = vj
    end
    return Star(c, V, P)
end

# ============================
# Star set getter functions
# ============================

"""
    center(X::Star)

Return the center of a star.

### Input

- `X` -- star

### Output

The center of the star.
"""
center(X::Star) = X.c

"""
    basis(X::Star)

Return the basis vectors of a star.

### Input

- `X` -- star

### Output

A matrix where each column is a basis vector of the star.
"""
basis(X::Star) = X.V

"""
    predicate(X::Star)

Return the predicate of a star.

### Input

- `X` -- star

### Output

A polyhedral set representing the predicate of the star.
"""
predicate(X::Star) = X.P

"""
    dim(X::Star)

Return the dimension of a star.

### Input

- `X` -- star set

### Output

The ambient dimension of the star set.
"""
dim(X::Star) = length(X.c)

# =====================================
# Support function and support vector
# =====================================

"""
    ρ(ℓ, X::Star)
Evaluate the support function of a star in a given direction.
### Input
- `ℓ` -- direction
- `X` -- star set
### Output
The support function in the given direction.
"""
function ρ(ℓ, X::Star)
    c = center(X)
    V = basis(X)
    P = predicate(X)
    return dot(ℓ, c) + ρ(transpose(V) * ℓ, P)
end
