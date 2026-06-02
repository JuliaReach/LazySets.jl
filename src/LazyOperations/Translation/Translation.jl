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

include("an_element.jl")
include("center.jl")
include("concretize.jl")
include("constraints_list.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("in.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
