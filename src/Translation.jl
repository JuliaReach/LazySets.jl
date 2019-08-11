import Base: isempty

export Translation,
       an_element,
       constraints_list,
       linear_map

"""
    Translation{N<:Real, VN<:AbstractVector{N}, S<:LazySet{N}} <: LazySet{N}

Type that represents a lazy translation.

The translation of set `X` along vector `v` is the map:

```math
x ↦ x + v,\\qquad x ∈ X
```

A translation is a special case of an affine map ``A x + b, x ∈ X`` where the
linear map ``A`` is the identity matrix and the translation vector ``b = v``.

### Fields

- `X` -- convex set
- `v` -- vector that defines the translation

### Example

```jldoctest translation
julia> X = BallInf([2.0, 2.0, 2.0], 1.0);

julia> v = [1.0, 0.0, 0.0]; # translation along dimension 1

julia> tr = Translation(X, v);

julia> typeof(tr)
Translation{Float64,Array{Float64,1},BallInf{Float64}}

julia> tr.X
BallInf{Float64}([2.0, 2.0, 2.0], 1.0)

julia> tr.v
3-element Array{Float64,1}:
 1.0
 0.0
 0.0
```
The sum operator `+` is overloaded to create translations:

```jldoctest translation
julia> X + v == Translation(X, v)
true
```
And so does the Minkowski sum operator, `⊕`:

```jldoctest translation
julia> X ⊕ v == Translation(X, v)
true
```

The translation of a translation is performed immediately:

```jldoctest translation
julia> tr = (X+v)+v
Translation{Float64,Array{Float64,1},BallInf{Float64}}(BallInf{Float64}([2.0, 2.0, 2.0], 1.0), [2.0, 0.0, 0.0])

julia> tr.v
3-element Array{Float64,1}:
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
`ρ` respectively:

```jldoctest translation
julia> σ([1.0, 0.0, 0.0], tr)
3-element Array{Float64,1}:
 5.0
 3.0
 3.0

julia> ρ([1.0, 0.0, 0.0], tr)
5.0
```
See the docstring of each of these functions for details.

The `an_element` function is useful to obtain an element of a translation:

```jldoctest translation
julia> e = an_element(tr)
3-element Array{Float64,1}:
 4.0
 2.0
 2.0
```

The lazy linear map of a translation is again a translation, since the following
simplification rule applies: ``M * (X⊕v) = (M*X) ⊕ (M*v)``:

```jldoctest translation
julia> using LinearAlgebra: I

julia> Q = Matrix(2.0I, 3, 3) * tr;

julia> Q isa Translation && Q.v == 2 * tr.v
true
```

Use the `isempty` method to query if the translation is empty; it falls back
to the `isempty` method of the wrapped set:

```jldoctest translation
julia> isempty(tr)
false
```

The list of constraints of the translation of a polyhedron (in general, a set
whose `constraints_list` is available) can be computed from a lazy translation:

```jldoctest translation
julia> constraints_list(tr)
6-element Array{HalfSpace{Float64,VN} where VN<:AbstractArray{Float64,1},1}:
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([1.0, 0.0, 0.0], 5.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([0.0, 1.0, 0.0], 3.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([0.0, 0.0, 1.0], 3.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([-1.0, 0.0, 0.0], -3.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([0.0, -1.0, 0.0], -1.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([0.0, 0.0, -1.0], -1.0)
```
"""
struct Translation{N<:Real, VN<:AbstractVector{N}, S<:LazySet{N}} <: LazySet{N}
    X::S
    v::VN

    # default constructor with dimension check
    function Translation(X::S, v::VN) where {N, VN<:AbstractVector{N}, S<:LazySet{N}}
        @assert dim(X) == length(v) "cannot create a translation of a set of dimension $(dim(X)) " *
                                    "along a vector of length $(length(v))"
        return new{N, VN, S}(X, v)
    end
end

# constructor from a Translation: perform the translation immediately
Translation(tr::Translation{N}, v::AbstractVector{N}) where {N<:Real} =
    Translation(tr.X, tr.v + v)

"""
    +(X::LazySet, v::AbstractVector)

Convenience constructor for a translation.

### Input

- `X` -- convex set
- `v` -- vector

### Output

The symbolic translation of ``X`` along vector ``v``.
"""
+(X::LazySet, v::AbstractVector) = Translation(X, v)

# translation from the left
+(v::AbstractVector, X::LazySet) = Translation(X, v)

"""
    ⊕(X::LazySet, v::AbstractVector)

Unicode alias constructor ⊕ (`oplus`) for the lazy translation operator.
"""
⊕(X::LazySet, v::AbstractVector) = Translation(X, v)

# translation from the left
⊕(v::AbstractVector, X::LazySet) = Translation(X, v)

# ============================
# LazySet interface functions
# ============================

"""
    dim(tr::Translation)::Int

Return the dimension of a translation.

### Input

- `tr` -- translation

### Output

The dimension of a translation.
"""
function dim(tr::Translation)::Int
    return length(tr.v)
end


"""
    σ(d::AbstractVector{N}, tr::Translation{N}) where {N<:Real}

Return the support vector of a translation.

### Input

- `d`  -- direction
- `tr` -- translation

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
function σ(d::AbstractVector{N}, tr::Translation{N}) where {N<:Real}
    return tr.v + σ(d, tr.X)
end

"""
    ρ(d::AbstractVector{N}, tr::Translation{N}) where {N<:Real}

Return the support function of a translation.

### Input

- `d`  -- direction
- `tr` -- translation

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector{N}, tr::Translation{N}) where {N<:Real}
    return dot(d, tr.v) + ρ(d, tr.X)
end

"""
    LinearMap(M::AbstractMatrix{N}, tr::Translation{N}) where {N<:Real}

Return the lazy linear map of a translation.

### Input

- `M`  -- matrix
- `tr` -- translation

### Output

The translation defined by the linear map.

### Notes

This method defines the simplification rule: ``M * (X⊕v) = (M*X) ⊕ (M*v)``,
returning a new translation.
"""
function LinearMap(M::AbstractMatrix{N}, tr::Translation{N}) where {N<:Real}
    return Translation(M * tr.X, M * tr.v)
end

"""
    an_element(tr::Translation)

Return some element of a translation.

### Input

- `tr` -- translation

### Output

An element in the translation.

### Notes

This function first asks for `an_element` function of the wrapped set, then
translates this element according to the given translation vector.
"""
function an_element(tr::Translation)
    return an_element(tr.X) + tr.v
end

"""
    isempty(tr::Translation)::Bool

Return if a translation is empty or not.

### Input

- `tr` -- translation

### Output

`true` iff the wrapped set is empty.
"""
function isempty(tr::Translation)::Bool
    return isempty(tr.X)
end

"""
    constraints_list(tr::Translation{N}, ::Val{true}) where {N<:Real}

Return the list of constraints of the translation of a set.

### Input

- `tr` -- lazy translation of a polyhedron

### Output

The list of constraints of the translation.

### Notes

We assume that the set wrapped by the lazy translation `X` offers a method
`constraints_list(⋅)`.

### Algorithm

Let the translation be defined by the set of points `y` such that `y = x + v` for
all `x ∈ X`. Then, each defining halfspace `a⋅x ≤ b` is transformed to
`a⋅y ≤ b + a⋅v`.
"""
function constraints_list(tr::Translation{N}, ::Val{true}) where {N<:Real}
    constraints_X = constraints_list(tr.X)
    constraints_TX = similar(constraints_X)
    @inbounds for (i, ci) in enumerate(constraints_X)
        constraints_TX[i] = HalfSpace(ci.a, ci.b + dot(ci.a, tr.v))
    end
    return constraints_TX
end

function constraints_list(tr::Translation{N}) where {N<:Real}
    has_constraints = applicable(constraints_list, tr.X)
    return constraints_list(tr, Val(has_constraints))
end

function constraints_list(tr::Translation{N}, ::Val{false}) where {N<:Real}
    throw(MethodError("this function requires that the `constraints_list` method is applicable"))
end

"""
    ∈(x::AbstractVector{N}, tr::Translation{N})::Bool where {N<:Real}

Check whether a given point is contained in the translation of a convex set.

### Input

- `x`  -- point/vector
- `tr` -- translation of a convex set

### Output

`true` iff ``x ∈ tr``.

### Algorithm

This implementation relies on the set membership function for the wrapped set
`tr.X`, since ``x ∈ X ⊕ v`` iff ``x - v ∈ X``.
"""
function ∈(x::AbstractVector{N}, tr::Translation{N})::Bool where {N<:Real}
    return x - tr.v ∈ tr.X
end

"""
    linear_map(M::AbstractMatrix{N}, tr::Translation{N}) where {N<:Real}

Concrete linear map of a polyhedron in constraint representation.

### Input

- `M`  -- matrix
- `tr` -- translation of a convex set

### Output

A concrete set corresponding to the linear map.
The type of the result depends on the type of the set wrapped by `tr`.

### Algorithm

We compute `translate(linear_map(M, tr.X), M * tr.v)`.
"""
function linear_map(M::AbstractMatrix{N}, tr::Translation{N}) where {N<:Real}
    @assert dim(tr) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
        "applied to a set of dimension $(dim(tr))"

    return translate(linear_map(M, tr.X), M * tr.v)
end
