import Base: *

"""
    LinearMap{N,S<:LazySet{N},NM,
              MAT<:Union{AbstractMatrix{NM},AbstractMatrixZonotope{NM}}} <: AbstractAffineMap{N,S}

Type that represents a linear transformation ``M⋅X`` of a set ``X``.

### Fields

- `M` -- linear map; can be a concrete matrix (`AbstractMatrix`) or a set-valued matrix (`AbstractMatrixZonotope`)
- `X` -- set

### Notes

This type is parametric in the elements of the linear map, `NM`, which is
independent of the numeric type of the wrapped set (`N`).
Typically `NM = N`, but there may be exceptions, e.g., if `NM` is an interval
that holds numbers of type `N`, where `N` is a floating point number type such
as `Float64`.

The linear map preserves convexity: if `X` is convex, then any linear map of `X`
is convex as well.

### Examples

For the examples we create a ``3×2`` matrix and a two-dimensional unit square.

```jldoctest constructors
julia> M = [1 2; 1 3; 1 4]; X = BallInf([0, 0], 1);
```

The function ``*`` can be used as an alias to construct a `LinearMap` object.

```jldoctest constructors
julia> lm = LinearMap(M, X)
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([1 2; 1 3; 1 4], BallInf{Int64, Vector{Int64}}([0, 0], 1))

julia> lm2 = M * X
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([1 2; 1 3; 1 4], BallInf{Int64, Vector{Int64}}([0, 0], 1))

julia> lm == lm2
true
```

For convenience, `M` does not need to be a matrix; we also allow to use vectors
(interpreted as an ``n×1`` matrix) and `UniformScaling`s resp. scalars
(interpreted as a scaling, i.e., a scaled identity matrix).
Scaling by ``1`` is ignored.

```jldoctest constructors
julia> using LinearAlgebra: I

julia> Y = BallInf([0], 1);  # one-dimensional interval

julia> [2, 3] * Y
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([2; 3;;], BallInf{Int64, Vector{Int64}}([0], 1))

julia> lm3 = 2 * X
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, SparseArrays.SparseMatrixCSC{Int64, Int64}}(sparse([1, 2], [1, 2], [2, 2], 2, 2), BallInf{Int64, Vector{Int64}}([0, 0], 1))

julia> 2I * X == lm3
true

julia> 1I * X == X
true
```

Applying a linear map to a `LinearMap` object combines the two maps into a
single `LinearMap` instance.
Again we can make use of the conversion for convenience.

```jldoctest constructors
julia> B = transpose(M); B * lm
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([3 9; 9 29], BallInf{Int64, Vector{Int64}}([0, 0], 1))

julia> B = [3, 4, 5]; B * lm
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([12 38], BallInf{Int64, Vector{Int64}}([0, 0], 1))

julia> B = 2; B * lm
LinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([2 4; 2 6; 2 8], BallInf{Int64, Vector{Int64}}([0, 0], 1))
```

The application of a `LinearMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest constructors
julia> M * ZeroSet{Int}(2)
ZeroSet{Int64}(3)

julia> M * EmptySet{Int}(2)
EmptySet{Int64}(3)
```
"""
struct LinearMap{N,S<:LazySet{N},NM,
                 MAT<:Union{AbstractMatrix{NM},AbstractMatrixZonotope{NM}}} <:
       AbstractAffineMap{N,S}
    M::MAT
    X::S

    # default constructor with dimension check
    function LinearMap(M::MAT,
                       X::S) where {N,S<:LazySet{N},NM,
                                    MAT<:Union{AbstractMatrix{NM},AbstractMatrixZonotope{NM}}}
        @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                     "be applied to a set of dimension $(dim(X))"
        return new{N,S,NM,MAT}(M, X)
    end
end

"""
```
    *(M::Union{AbstractMatrix, UniformScaling, AbstractVector, Real, AbstractMatrixZonotope},
      X::LazySet)
```

Alias to create a `LinearMap` object.

### Input

- `M` -- matrix or matrix zonotope
- `X` -- set

### Output

A lazy linear map, i.e., a `LinearMap` instance.
"""
function *(M::Union{AbstractMatrix,UniformScaling,AbstractVector,Real,AbstractMatrixZonotope},
           X::LazySet)
    return LinearMap(M, X)
end

# scaling from the right
function *(X::LazySet, M::Real)
    return LinearMap(M, X)
end

# convenience constructor from a vector
function LinearMap(v::AbstractVector, X::LazySet)
    return _LinearMap_vector(v, X)
end

function _LinearMap_vector(v, X)
    n = dim(X)
    m = length(v)
    if n == m
        M = reshape(v, 1, length(v))
    else
        M = reshape(v, length(v), 1)
    end
    return LinearMap(M, X)
end

# convenience constructor from a UniformScaling
function LinearMap(M::UniformScaling, X::LazySet)
    if isone(M.λ)
        return X
    end
    return LinearMap(Diagonal(fill(M.λ, dim(X))), X)
end

# convenience constructor from a scalar
function LinearMap(α::Real, X::LazySet)
    n = dim(X)
    return LinearMap(sparse(α * I, n, n), X)
end

# combine two linear maps into a single linear map
function LinearMap(M::AbstractMatrix, lm::LinearMap)
    return LinearMap(M * lm.M, lm.X)
end

# disambiguation
function LinearMap(v::AbstractVector, lm::LinearMap)
    return _LinearMap_vector(v, lm)
end

# more efficient versions when combining `LinearMap`s
function LinearMap(M::UniformScaling, lm::LinearMap)
    if isone(M.λ)
        return lm
    end
    return LinearMap(M.λ * lm.M, lm.X)
end

function LinearMap(α::Real, lm::LinearMap)
    return LinearMap(α * lm.M, lm.X)
end

# ZeroSet is "almost absorbing" for LinearMap (only the dimension changes)
function LinearMap(M::AbstractMatrix, Z::ZeroSet)
    N = promote_type(eltype(M), eltype(Z))
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                 "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

# EmptySet is "almost absorbing" for LinearMap (only the dimension changes)
function LinearMap(M::AbstractMatrix, ∅::EmptySet)
    N = promote_type(eltype(M), eltype(∅))
    @assert dim(∅) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                 "be applied to a set of dimension $(dim(∅))"
    return EmptySet{N}(size(M, 1))
end

function matrix(lm::LinearMap)
    return lm.M
end

function vector(lm::LinearMap{N}) where {N}
    return spzeros(N, dim(lm))
end

function set(lm::LinearMap)
    return lm.X
end

"""
    Projection(X::LazySet{N}, variables::AbstractVector{Int}) where {N}

Return a lazy projection of a set.

### Input

- `X`         -- set
- `variables` -- variables of interest

### Output

A lazy `LinearMap` that corresponds to projecting `X` along the given variables
`variables`.

### Examples

The projection of a three-dimensional cube into the first two coordinates:

```jldoctest Projection
julia> B = BallInf([1.0, 2, 3], 1.0)
BallInf{Float64, Vector{Float64}}([1.0, 2.0, 3.0], 1.0)

julia> Bproj = Projection(B, [1, 2])
LinearMap{Float64, BallInf{Float64, Vector{Float64}}, Float64, SparseArrays.SparseMatrixCSC{Float64, Int64}}(sparse([1, 2], [1, 2], [1.0, 1.0], 2, 3), BallInf{Float64, Vector{Float64}}([1.0, 2.0, 3.0], 1.0))

julia> isequivalent(Bproj, BallInf([1.0, 2], 1.0))
true
```
"""
function Projection(X::LazySet{N}, variables::AbstractVector{Int}) where {N}
    M = projection_matrix(variables, dim(X), N)
    return LinearMap(M, X)
end

include("an_element.jl")
include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isoperationtype.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
