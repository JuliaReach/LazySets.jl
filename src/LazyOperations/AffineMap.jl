export AffineMap

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
AffineMap{Int64,BallInf{Int64,Array{Int64,1}},Int64,Array{Int64,2},Array{Int64,1}}([1 2; 1 3; 1 4], BallInf{Int64,Array{Int64,1}}([0, 0], 1), [1, 2, 3])
```

For convenience, `A` does not need to be a matrix but we also allow to use
`UniformScaling`s resp. scalars (interpreted as a scaling, i.e., a scaled
identity matrix).
Scaling by ``1`` is ignored and simplified to a pure `Translation`.

```jldoctest constructors
julia> using LinearAlgebra

julia> am = AffineMap(2I, X, b2)
AffineMap{Int64,BallInf{Int64,Array{Int64,1}},Int64,Diagonal{Int64,Array{Int64,1}},Array{Int64,1}}([2 0; 0 2], BallInf{Int64,Array{Int64,1}}([0, 0], 1), [1, 2])

julia> AffineMap(2, X, b2) == am
true

julia> AffineMap(1, X, b2)
Translation{Int64,Array{Int64,1},BallInf{Int64,Array{Int64,1}}}(BallInf{Int64,Array{Int64,1}}([0, 0], 1), [1, 2])
```

Applying a linear map to an `AffineMap` object combines the two maps into a new
`AffineMap` instance.
Again we can make use of the conversion for convenience.

```jldoctest constructors
julia> B = [2 0; 0 2]; am2 = B * am
AffineMap{Int64,BallInf{Int64,Array{Int64,1}},Int64,Array{Int64,2},Array{Int64,1}}([4 0; 0 4], BallInf{Int64,Array{Int64,1}}([0, 0], 1), [2, 4])

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


function matrix(am::AffineMap)
    return am.M
end

function vector(am::AffineMap{N}) where {N<:Real}
    return am.v
end

function set(am::AffineMap)
    return am.X
end
