export Star,
       basis,
       predicate,
       intersection!

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

isoperationtype(::Type{<:Star}) = false

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

- `X` -- star

### Output

The ambient dimension of a star.
"""
function dim(X::Star)
    return length(X.c)
end

"""
    σ(d::AbstractVector, X::Star)

Return a support vector of a star.

### Input

- `d` -- direction
- `X` -- star

### Output

A support vector in the given direction.
"""
function σ(d::AbstractVector, X::Star)
    A = basis(X)
    return A * σ(At_mul_B(A, d), predicate(X)) + center(X)
end

"""
    ρ(d::AbstractVector, X::Star)

Evaluate the support function of a star.

### Input

- `d` -- direction
- `X` -- star

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector, X::Star)
    return ρ(At_mul_B(basis(X), d), predicate(X)) + dot(d, center(X))
end

"""
    an_element(X::Star)

Return some element of a star.

### Input

- `X` -- star

### Output

An element of the star.

### Algorithm

We apply the affine map to the result of `an_element` on the predicate.
"""
function an_element(X::Star)
    return basis(X) * an_element(predicate(X)) + center(X)
end

"""
    isempty(X::Star)

Check whether a star is empty.

### Input

- `X` -- star

### Output

`true` iff the predicate is empty.
"""
function isempty(X::Star)
    return isempty(predicate(X))
end

"""
    isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)

Check whether a star is bounded.

### Input

- `X`        -- star
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the basis matrix is invertible)

### Output

`true` iff the star is bounded.

### Algorithm

See [`isbounded(::AbstractAffineMap)`](@ref).
"""
function isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)
    am = convert(STAR, X)
    return isbounded(am; cond_tol=cond_tol)
end

"""
    ∈(v::AbstractVector, X::Star)

Check whether a given point is contained in a star.

### Input

- `v` -- point/vector
- `X` -- star

### Output

`true` iff ``v ∈ X``.

### Algorithm

The implementation is identical to
[`∈(::AbstractVector, ::AbstractAffineMap)`](@ref).
"""
function ∈(x::AbstractVector, X::Star)
    return basis(X) \ (x - center(X)) ∈ predicate(X)
end

"""
    vertices_list(X::Star; apply_convex_hull::Bool=true)

Return the list of vertices of a star.

### Input

- `X`                 -- star
- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices of the star

### Output

A list of vertices.

### Algorithm

See [`vertices_list(::AbstractAffineMap)`](@ref).
"""
function vertices_list(X::Star; apply_convex_hull::Bool=true)
    am = convert(STAR, X)
    return vertices_list(am; apply_convex_hull=apply_convex_hull)
end

"""
    constraints_list(X::Star)

Return the list of constraints of a star.

### Input

- `X` -- star

### Output

The list of constraints of the star.

### Algorithm

See [`constraints_list(::AbstractAffineMap)`](@ref).
"""
function constraints_list(X::Star)
    am = convert(STAR, X)
    return constraints_list(am)
end

"""
    linear_map(M::AbstractMatrix, X::Star)

Return the linear map of a star.

### Input

- `M` -- matrix
- `X` -- star

### Output

The star obtained by applying `M` to `X`.
"""
function linear_map(M::AbstractMatrix, X::Star)
    c′ = M * X.c
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end

"""
    affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)

Return the affine map of a star.

### Input

- `M` -- matrix
- `X` -- star
- `v` -- vector

### Output

The star obtained by applying the affine map with matrix `M` and displacement
`v` to `X`.
"""
function affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)
    c′ = M * X.c + v
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end

"""
    rand(::Type{Star}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random star.

### Input

- `Star` -- type for dispatch
- `N`    -- (optional, default: `Float64`) numeric type
- `dim`  -- (optional, default: 2) dimension
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding
- `P`    -- (optional, default: a random `HPolytope`) predicate

### Output

A random star.

### Algorithm

By default we generate a random `HPolytope` of dimension `dim` as predicate.
Alternatively the predicate can be passed.

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Star};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              P::AbstractPolyhedron=rand(HPolytope; N=N, dim=dim, rng=rng, seed=seed))
    rng = reseed!(rng, seed)
    c = randn(rng, N, dim)
    V = randn(rng, N, dim, LazySets.dim(P))
    return Star(c, V, P)
end
