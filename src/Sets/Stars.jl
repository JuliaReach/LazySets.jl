export center,
       basis,
       predicate

"""
    AbstractStarSet{N} <: LazySet{N}

Abstract supertype for all star set types.
"""
abstract type AbstractStarSet{N} <: LazySet{N} end

isoperationtype(::Type{<:AbstractStarSet}) = false
isconvextype(::Type{<:AbstractStarSet}) = false

function linear_map(M::AbstractMatrix, X::AbstractStarSet)
    # .......
end

function affine_map(M::AbstractMatrix, X::AbstractStarSet, v::AbstractVector)
    # .......
end

function intersection(X::AbstractStarSet, H::HalfSpace)
    # .......
end

# symmetric method
intersection(H::HalfSpacev, X::AbstractStarSet) = intersection(X, H)

"""
    GeneralizedStar{N, VN<:AbstractVector{N}, FT} <: LazySet{N}

Type that represents a generalized star set.

### Fields

- `c` -- vector that represents the center
- `V` -- list of vectors that represents the basis
- `P` -- function that represents the predicate

### Notes

This type does not make any a-priori assumption on the predicate's nature.
On the contrary, the type `Star` assumes that the predicate is polyhedral.

For definitions and an application of star sets in the context of verification
of neural networks see [1].

### References

- [1] Duggirala, Parasara Sridhar, and Mahesh Viswanathan.
      *Parsimonious, simulation based verification of linear systems.*
      International Conference on Computer Aided Verification. Springer, Cham, 2016.
"""
struct GeneralizedStar{N, VN<:AbstractVector{N}, FT} <: AbstractStarSet{N}
    c::VN # center
    V::Vector{VN} # basis
    P::FT # predicate
end

# getter functions
center(X::GeneralizedStar) = X.s
basis(X::GeneralizedStar) = X.V
predicate(X::GeneralizedStar) = X.P

"""
    Star{N, VN<:AbstractVector{N}} <: AbstractStarSet{N}

Type that represents a star set, i.e. a generalized star with a polyhedral predicate.

### Fields

- `c` -- vector that represents the center
- `V` -- list of vectors that represents the basis
- `P` -- polyhedron that represents the predicate

### Notes

In particular, the predicate function is implemented as a conjunction of linear
constraints. This definition can be found in reference [1].

### Examples

This example is drawn from Example 1 in [1]. Consider the two-dimensional plane
``\\mathbb{R}^2``. Let

```jldoctest star_constructor
julia> V = [[1.0, 0.0], [0.0, 1.0]]

```
be the set of unit vectors along the two axes, and take

```jldoctest star_constructor
julia> c = [3.0, 3.0]
```
be the center of the star set. Then, the star set `X = Star(c, V, P)` defines
the set:

```jldoctest star_constructor
julia>
```

### References

- [1] Duggirala, Parasara Sridhar, and Mahesh Viswanathan.
      *Parsimonious, simulation based verification of linear systems.*
      International Conference on Computer Aided Verification. Springer, Cham, 2016.
"""
struct Star{N, VN<:AbstractVector{N}} <: AbstractStarSet{N}
    c::VN # center
    V::Vector{VN} # basis
    P::HPolyhedron{N, VN} # predicate
end

# getter functions
center(X::Star) = X.s
basis(X::Star) = X.V
predicate(X::Star) = X.P
