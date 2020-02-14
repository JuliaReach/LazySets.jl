import Base: isempty

export AffineMap,
       an_element,
       isempty,
       isbounded,
       ∈,
       vertices_list,
       constraints_list,
       linear_map

"""
    AffineMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM},
              VN<:AbstractVector{NM}} <: AbstractAffineMap{N, S}

Type that represents an affine transformation ``M⋅X ⊕ v`` of a convex set ``X``.

### Fields

- `M` -- matrix/linear map
- `X` -- convex set
- `v` -- translation vector

### Notes

An affine map is the composition of a linear map and a translation. This type is
parametric in the coefficients of the linear map, `NM`, which may be different from
the numeric type of the wrapped set (`N`). However, the numeric type of the
translation vector should be `NM`.

### Examples

For the examples we create a ``3×2`` matrix, a two-dimensional unit square, and
a three-dimensional vector.
Then we combine them in an `AffineMap`.

```jldoctest constructors
julia> A = [1 2; 1 3; 1 4]; X = BallInf([0, 0], 1); b2 = [1, 2]; b3 = [1, 2, 3];

julia> AffineMap(A, X, b3)
AffineMap{Int64,BallInf{Int64},Int64,Array{Int64,2},Array{Int64,1}}([1 2; 1 3; 1 4], BallInf{Int64}([0, 0], 1), [1, 2, 3])
```

For convenience, `A` does not need to be a matrix but we also allow to use
`UniformScaling`s resp. scalars (interpreted as a scaling, i.e., a scaled
identity matrix).
Scaling by ``1`` is ignored and simplified to a pure `Translation`.

```jldoctest constructors
julia> using LinearAlgebra

julia> am = AffineMap(2I, X, b2)
AffineMap{Int64,BallInf{Int64},Int64,Diagonal{Int64,Array{Int64,1}},Array{Int64,1}}([2 0; 0 2], BallInf{Int64}([0, 0], 1), [1, 2])

julia> AffineMap(2, X, b2) == am
true

julia> AffineMap(1, X, b2)
Translation{Int64,Array{Int64,1},BallInf{Int64}}(BallInf{Int64}([0, 0], 1), [1, 2])
```

Applying a linear map to an `AffineMap` object combines the two maps into a new
`AffineMap` instance.
Again we can make use of the conversion for convenience.

```jldoctest constructors
julia> B = [2 0; 0 2]; am2 = B * am
AffineMap{Int64,BallInf{Int64},Int64,Array{Int64,2},Array{Int64,1}}([4 0; 0 4], BallInf{Int64}([0, 0], 1), [2, 4])

julia> 2 * am == am2
true
```

The application of an `AffineMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest constructors
julia> AffineMap(A, ZeroSet{Int}(2), b3)
Singleton{Int64,Array{Int64,1}}([1, 2, 3])

julia> AffineMap(A, EmptySet{Int}(2), b3)
EmptySet{Int64}(2)
```
"""
struct AffineMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM},
                 VN<:AbstractVector{NM}} <: AbstractAffineMap{N, S}
    M::MAT
    X::S
    v::VN

    # default constructor with dimension match check
    function AffineMap(M::MAT, X::S, v::VN) where {N<:Real, S<:LazySet{N}, NM,
                                                   MAT<:AbstractMatrix{NM},
                                                   VN<:AbstractVector{NM}}

        @assert dim(X) == size(M, 2) "a matrix of size $(size(M)) cannot be " *
            "applied to a set of dimension $(dim(X))"

        @assert size(M, 1) == length(v) "a map with output dimension " *
            "$(size(M, 1)) is incompatible with a translation vector of " *
            "dimension $(length(v))"

        return new{N, S, NM, MAT, VN}(M, X, v)
    end
end

isoperationtype(::Type{<:AffineMap}) = true
isconvextype(::Type{<:AffineMap{N, S}}) where {N, S} = isconvextype(S)

# convenience constructor from a UniformScaling
function AffineMap(M::UniformScaling, X::LazySet, v::AbstractVector)
    return AffineMap(M.λ, X, v)
end

# convenience constructor from a scalar
function AffineMap(α::N, X::LazySet, v::AbstractVector) where {N<:Real}
    if α == one(N)
        return Translation(X, v)
    end
    return AffineMap(Diagonal(fill(α, length(v))), X, v)
end

# simplification for a LinearMap
for MAP in (:AbstractMatrix, :Real)
    @eval begin
        function LinearMap(map::$MAP, am::AffineMap)
            return AffineMap(map * am.M, am.X, map * am.v)
        end
    end
end

# ZeroSet is "almost absorbing" for the linear map (only the dimension changes)
# such that only the translation vector remains
function AffineMap(M::AbstractMatrix{N}, Z::ZeroSet{N}, v::AbstractVector{N}
                  ) where {N<:Real}
    @assert dim(Z) == size(M, 2) "a matrix of size $(size(M)) cannot be " *
        "applied to a set of dimension $(dim(Z))"

    @assert size(M, 1) == length(v) "a map with output dimension " *
        "$(size(M, 1)) is incompatible with a translation vector of " *
        "dimension $(length(v))"

    return Singleton(v)
end

# EmptySet is absorbing for AffineMap
function AffineMap(M::AbstractMatrix{N}, ∅::EmptySet{N}, v::AbstractVector{N}
                  ) where {N<:Real}
    return ∅
end


# --- AbstractAffineMap interface functions ---


function get_A(am::AffineMap)
    return am.M
end

function get_b(am::AffineMap{N}) where {N<:Real}
    return am.v
end

function get_X(am::AffineMap)
    return am.X
end

# ============================
# LazySet interface functions
# ============================

"""
    dim(am::AffineMap)

Return the dimension of an affine map.

### Input

- `am` -- affine map

### Output

The dimension of an affine map.
"""
function dim(am::AffineMap)
    return length(am.v)
end

"""
    σ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}

Return the support vector of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}
    return am.M * σ(_At_mul_B(am.M, d), am.X) + am.v
end

"""
    ρ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}

Return the support function of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}
    return ρ(_At_mul_B(am.M, d), am.X) + dot(d, am.v)
end

"""
    an_element(am::AffineMap)

Return some element of an affine map.

### Input

- `am` -- affine map

### Output

An element of the affine map. It relies on the `an_element` function of the
wrapped set.
"""
function an_element(am::AffineMap)
    return am.M * an_element(am.X) + am.v
end

"""
    isempty(am::AffineMap)

Return whether an affine map is empty or not.

### Input

- `am` -- affine map

### Output

`true` iff the wrapped set is empty and the affine vector is empty.
"""
function isempty(am::AffineMap)
    return isempty(am.X)
end

"""
    isbounded(am::AffineMap; cond_tol::Number=DEFAULT_COND_TOL)

Determine whether an affine map is bounded.

### Input

- `am`       -- affine map
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)

### Output

`true` iff the affine map is bounded.

### Algorithm

We first check if the matrix is zero or the wrapped set is bounded.
If not, we perform a sufficient check whether the matrix is invertible.
If the matrix is invertible, then the map being bounded is equivalent to the
wrapped set being bounded, and hence the map is unbounded.
Otherwise, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(am::AffineMap; cond_tol::Number=DEFAULT_COND_TOL)
    if iszero(am.M) || isbounded(am.X)
        return true
    end
    if isinvertible(am.M; cond_tol=cond_tol)
        return false
    end
    return isbounded_unit_dimensions(am)
end

"""
    ∈(x::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}

Check whether a given point is contained in the affine map of a convex set.

### Input

- `x`  -- point/vector
- `am` -- affine map of a convex set

### Output

`true` iff ``x ∈ am``.

### Algorithm

Note that ``x ∈ M⋅S ⊕ v`` iff ``M^{-1}⋅(x - v) ∈ S``.
This implementation does not explicitly invert the matrix, which is why it also
works for non-square matrices.

### Examples

```jldoctest
julia> am = AffineMap([2.0 0.0; 0.0 1.0], BallInf([1., 1.], 1.), [-1.0, -1.0]);

julia> [5.0, 1.0] ∈ am
false

julia> [3.0, 1.0] ∈ am
true
```

An example with a non-square matrix:

```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> [0.5, 0.5] ∈ M*B
true
```
"""
function ∈(x::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}
    return ∈(am.M \ (x - am.v), am.X)
end

"""
    vertices_list(am::AffineMap{N}; [apply_convex_hull]::Bool) where {N<:Real}

Return the list of vertices of a (polyhedral) affine map.

### Input

- `am`                -- affine map
- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices transformed by the
                         affine map

### Output

A list of vertices.

### Algorithm

This implementation computes all vertices of `X`, then transforms them through
the affine map, i.e. `x ↦ M*x + v` for each vertex `x` of `X`. By default, the
convex hull operation is taken before returning this list. For dimensions three
or higher, this operation relies on the functionality through the concrete
polyhedra library `Polyhedra.jl`.

If you are not interested in taking the convex hull of the resulting vertices under
the affine map, pass `apply_convex_hull=false` as a keyword argument.

Note that we assume that the underlying set `X` is polyhedral, either concretely
or lazily, i.e. there the function `vertices_list` should be applicable.
"""
function vertices_list(am::AffineMap{N};
                       apply_convex_hull::Bool=true) where {N<:Real}

    # collect vertices list of the wrapped set
    vlist_X = vertices_list(am.X)

    # create resulting vertices list
    vlist = [am.M * x + am.v for x in vlist_X]

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

"""
    constraints_list(am::AffineMap{N}) where {N<:Real}

Return the list of constraints of a (polyhedral) affine map.

### Input

- `am` -- affine map

### Output

The list of constraints of the affine map.

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

Falls back to the list of constraints of the translation of a lazy linear map.
"""
function constraints_list(am::AffineMap{N}) where {N<:Real}
    return _constraints_list_translation(LinearMap(am.M, am.X), am.v)
end

"""
    linear_map(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}

Return the linear map of a lazy affine map.

### Input

- `M`  -- matrix
- `am` -- affine map

### Output

A set corresponding to the linear map of the lazy affine map of a set.
"""
function linear_map(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}
     return translate(linear_map(M * am.M, am.X), M * am.v)
end
