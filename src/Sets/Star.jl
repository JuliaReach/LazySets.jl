export Star,
       center,
       basis,
       predicate

"""
    Star{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}} <: AbstractStar{N}

Type that represents a generalized star set where the predicate is polyhedral.

### Fields

- `c` -- vector that represents the center
- `V` -- matrix where each column corresponds to a basis vector
- `P` -- polyhedral set that represents the predicate

### Notes

The predicate function is implemented as a conjunction of linear constraints.
This definition can be found in reference [1].

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
struct Star{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, PT<:AbstractPolyhedron{N}} <: AbstractStar{N}
    c::VN # center
    V::MN # basis
    P::PT # predicate
    # TODO INNER CONSTRUCTOR WITH DIMENSION CHECK
end

# getter functions
center(X::Star) = X.c
basis(X::Star) = X.V
predicate(X::Star) = X.P

dim(X::Star) = length(X.c)

function linear_map!(M::AbstractMatrix, X::Star, c′, V′)
    @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot be applied " *
                                 "to a star of dimension $(dim(X))"

    c = center(X)
    V = basis(X)
    mul!(c′, M, c)
    mul!(V′, M, V)
    P = predicate(X)
    return Star(c′, V′, P)
end

function linear_map(M::AbstractMatrix, X::Star)
    @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot be applied " *
                                 "to a star of dimension $(dim(X))"

    c = center(X)
    V = basis(X)
    c′ = M * c
    V′ = M * V
    P = predicate(X)
    return Star(c′, V′, P)
end

# compute the affine map M * X + c′ overriding c′
function affine_map!(M::AbstractMatrix{N}, X::Star, c′, V′) where {N}
    @assert dim(X) == size(M, 2) == length(c′) "an affine map of size $(size(M)) and " *
    "translation of length $(length(c′)) cannot be applied to a star of dimension $(dim(X))"

    c = center(X)
    V = basis(X)
    mul!(c′, M, c, one(N), one(N))
    mul!(V′, M, V)
    P′ = predicate(X)
    return Star(c′, V′, P′)
end

function affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)
    @assert dim(X) == size(M, 2) == length(v) "an affine map of size $(size(M)) and " *
    "of dimension $(length(v)) cannot be applied to a star of dimension $(dim(X))"

    c = center(X)
    V = basis(X)
    P = predicate(X)
    c′ = M * c + v
    V′ = M * V
    return Star(c′, V′, P)
end

function addconstraint!(X::Star, H::HalfSpace)
    addconstraint!(X.P, H)
end

function remove_redundant_constraints(X::Star{N}; backend=default_lp_solver(N)) where {N}
    return remove_redundant_constraints(X.P, backend=backend)
end

function remove_redundant_constraints!(X::Star{N}; backend=default_lp_solver(N)) where {N}
    return remove_redundant_constraints!(X.P, backend=backend)
end

function _intersection!(c, V, P::AbstractPolyhedron, H::HalfSpace)
    a′ = transpose(V) * H.a
    b′ = H.b - dot(H.a, c)
    H′ = HalfSpace(a′, b′)
    addconstraint!(P, H′)
end

function intersection!(X::Star, H::HalfSpace)
    _intersection!(center(X), basis(X), predicate(X), H)
    return X
end

function intersection(X::Star, H::HalfSpace)
    intersection!(copy(X), H)
end

# symmetric methods
intersection(H::HalfSpace, X::Star) = intersection(X, H)
intersection!(H::HalfSpace, X::Star) = intersection!(X, H)
