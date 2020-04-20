export center,
       basis,
       predicate

"""
    Star{N, VN<:AbstractVector{N}} <: AbstractStarSet{N}

Type that represents a generalized star set.

### Fields

- `c` -- vector that represents the center
- `V` -- list of vectors that represents the basis
- `P` -- polyhedron that represents the predicate

### Notes

In particular, the predicate function is implemented as a conjunction of linear
constraints.

The notion defined here is taken from the *generalized star set representation*
in [1].

References:

- [1] Duggirala, Parasara Sridhar, and Mahesh Viswanathan.
      *Parsimonious, simulation based verification of linear systems.*
      International Conference on Computer Aided Verification. Springer, Cham, 2016.

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
Then, the generalized star set `Θ = Star(c, V, P)` defines the set

```jldoctest star_constructor
julia>
```
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
