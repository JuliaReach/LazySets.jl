export Translation

"""
    Translation{N, S<:LazySet{N}, VN<:AbstractVector{N}}
        <: AbstractAffineMap{N, S}

Type that represents a lazy translation.

The translation of set `X` along vector `v` is the map:

```math
x ↦ x + v,\\qquad x ∈ X
```

A translation is a special case of an affine map ``A x + b, x ∈ X`` where the
linear map ``A`` is the identity matrix and the translation vector ``b`` is
``v``.

### Fields

- `X` -- set
- `v` -- vector that defines the translation

### Notes

Translation preserves convexity: if `X` is convex, then any translation of `X`
is convex as well.

The convenience aliases `⊕` and `+` are also available. `⊕` can be typed by
`\\oplus<tab>`.

### Example

```jldoctest translation
julia> X = BallInf([2.0, 2.0, 2.0], 1.0);

julia> v = [1.0, 0.0, 0.0]; # translation along dimension 1

julia> tr = Translation(X, v);

julia> typeof(tr)
Translation{Float64, BallInf{Float64, Vector{Float64}}, Vector{Float64}}

julia> tr.X
BallInf{Float64, Vector{Float64}}([2.0, 2.0, 2.0], 1.0)

julia> tr.v
3-element Vector{Float64}:
 1.0
 0.0
 0.0
```

Both the sum operator `+` and the Minkowski-sum operator `⊕` are overloaded to
create translations:

```jldoctest translation
julia> X + v == X ⊕ v == Translation(X, v)
true
```

The translation of a translation is performed immediately:

```jldoctest translation
julia> tr = (X + v) + v
Translation{Float64, BallInf{Float64, Vector{Float64}}, Vector{Float64}}(BallInf{Float64, Vector{Float64}}([2.0, 2.0, 2.0], 1.0), [2.0, 0.0, 0.0])

julia> tr.v
3-element Vector{Float64}:
 2.0
 0.0
 0.0
```

The dimension of a translation is obtained with the `dim` function:

```jldoctest translation
julia> dim(tr)
3
```

For the support vector (resp. support function) along vector `d`, use `σ` and
`ρ`, respectively:

```jldoctest translation
julia> σ([1.0, 0.0, 0.0], tr)
3-element Vector{Float64}:
 5.0
 2.0
 2.0

julia> ρ([1.0, 0.0, 0.0], tr)
5.0
```

See the docstring of each of these functions for details.

The `an_element` function is useful to obtain an element of a translation:

```jldoctest translation
julia> e = an_element(tr)
3-element Vector{Float64}:
 4.0
 2.0
 2.0
```

The lazy linear map of a translation is an affine map, since the following
simplification rule applies: ``M * (X ⊕ v) = (M * X) ⊕ (M * v)``:

```jldoctest translation
julia> using LinearAlgebra: I

julia> M = Matrix(2.0I, 3, 3);

julia> Q = M * tr;

julia> Q isa AffineMap && Q.M == M && Q.X == tr.X && Q.v == 2 * tr.v
true
```

Use the `isempty` method to check whether the translation is empty:

```jldoctest translation
julia> isempty(tr)
false
```

The list of constraints of the translation of a polyhedral set (a set whose
`constraints_list` is available) can be computed from a lazy translation:

```jldoctest translation
julia> constraints_list(tr)
6-element Vector{HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}}:
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([1.0, 0.0, 0.0], 5.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([0.0, 1.0, 0.0], 3.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([0.0, 0.0, 1.0], 3.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([-1.0, 0.0, 0.0], -3.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([0.0, -1.0, 0.0], -1.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([0.0, 0.0, -1.0], -1.0)
```
"""
struct Translation{N,S<:LazySet{N},VN<:AbstractVector{N}} <: AbstractAffineMap{N,S}
    X::S
    v::VN

    # default constructor with dimension check
    function Translation(X::S, v::VN) where {N,VN<:AbstractVector{N},S<:LazySet{N}}
        @assert dim(X) == length(v) "cannot create a translation of a set of " *
                                    "dimension $(dim(X)) along a vector of length $(length(v))"

        return new{N,S,VN}(X, v)
    end
end

isoperationtype(::Type{<:Translation}) = true

isconvextype(::Type{Translation{N,S,VN}}) where {N,S,VN} = isconvextype(S)

# constructor from a Translation: perform the translation immediately
Translation(tr::Translation{N}, v::AbstractVector{N}) where {N} = Translation(tr.X, tr.v + v)

# the translation of a lazy linear map is a (lazy) affine map
Translation(lm::LinearMap, v::AbstractVector) = AffineMap(lm.M, lm.X, v)

# the linear map of a translation is a (lazy) affine map:
# M * (X ⊕ v) = (M * X) ⊕ (M * v)
LinearMap(M::AbstractMatrix, tr::Translation) = AffineMap(M, tr.X, M * tr.v)

# EmptySet is absorbing for Translation
Translation(∅::EmptySet, v::AbstractVector) = ∅

# Universe is absorbing for Translation
Translation(U::Universe, v::AbstractVector) = U

@commutative +(X::LazySet, v::AbstractVector) = Translation(X, v)

@commutative ⊕(X::LazySet, v::AbstractVector) = Translation(X, v)

function matrix(tr::Translation{N}) where {N}
    return Diagonal(ones(N, dim(tr)))
end

function vector(tr::Translation)
    return tr.v
end

function set(tr::Translation)
    return tr.X
end

"""
    σ(d::AbstractVector, tr::Translation)

Return a support vector of a translation.

### Input

- `d`  -- direction
- `tr` -- translation of a set

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
function σ(d::AbstractVector, tr::Translation)
    return tr.v + σ(d, tr.X)
end

"""
    ρ(d::AbstractVector, tr::Translation)

Evaluate the support function of a translation.

### Input

- `d`  -- direction
- `tr` -- translation of a set

### Output

The evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, tr::Translation)
    return dot(d, tr.v) + ρ(d, tr.X)
end

"""
    an_element(tr::Translation)

Return some element of a translation.

### Input

- `tr` -- translation of a set

### Output

An element in the translation.

### Notes

This function first asks for `an_element` of the wrapped set, then translates
this element according to the given translation vector.
"""
function an_element(tr::Translation)
    return an_element(tr.X) + tr.v
end

function isboundedtype(::Type{<:Translation{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
    constraints_list(tr::Translation)

Return a list of constraints of the translation of a set.

### Input

- `tr` -- translation of a polyhedron

### Output

A list of constraints of the translation.

### Notes

We assume that the set wrapped by the lazy translation `X` offers a method
`constraints_list(⋅)`.

### Algorithm

Let the translation be defined by the set of points `y` such that `y = x + v` for
all `x ∈ X`. Then, each defining halfspace `a⋅x ≤ b` is transformed to
`a⋅y ≤ b + a⋅v`.
"""
function constraints_list(tr::Translation)
    return _constraints_list_translation(tr.X, tr.v)
end

function _constraints_list_translation(X::LazySet, v::AbstractVector)
    constraints_X = constraints_list(X)
    constraints_TX = similar(constraints_X)
    @inbounds for (i, ci) in enumerate(constraints_X)
        constraints_TX[i] = HalfSpace(ci.a, ci.b + dot(ci.a, v))
    end
    return constraints_TX
end

"""
    ∈(x::AbstractVector, tr::Translation)

Check whether a given point is contained in the translation of a set.

### Input

- `x`  -- point/vector
- `tr` -- translation of a set

### Output

`true` iff ``x ∈ tr``.

### Algorithm

This implementation relies on the set-membership function for the wrapped set
`tr.X`, since ``x ∈ X ⊕ v`` iff ``x - v ∈ X``.
"""
function ∈(x::AbstractVector, tr::Translation)
    return x - tr.v ∈ tr.X
end

"""
    linear_map(M::AbstractMatrix, tr::Translation)

Concrete linear map of a translation.

### Input

- `M`  -- matrix
- `tr` -- translation of a set

### Output

A concrete set corresponding to the linear map.
The type of the result depends on the type of the set wrapped by `tr`.

### Algorithm

We compute `affine_map(M, tr.X, M * tr.v)`.
"""
function linear_map(M::AbstractMatrix, tr::Translation)
    @assert dim(tr) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                  "applied to a set of dimension $(dim(tr))"

    return affine_map(M, tr.X, M * tr.v)
end

function concretize(tr::Translation)
    return translate(concretize(tr.X), tr.v)
end

"""
    center(tr::Translation)

Return the center of the translation of a centrally-symmetric set.

### Input

- `tr` -- translation of a centrally-symmetric set

### Output

The translation of the center of the wrapped set by the translation vector.
"""
function center(tr::Translation)
    return center(tr.X) + tr.v
end

function translate(tr::Translation, x::AbstractVector)
    return Translation(translate(tr.X, x))
end

function isuniversal(tr::Translation)
    return isuniversal(tr.X)
end
