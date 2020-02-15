import Base: *, ∈, isempty

export LinearMap,
       an_element,
       constraints_list

"""
    LinearMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} <: LazySet{N}

Type that represents a linear transformation ``M⋅S`` of a convex set ``S``.

### Fields

- `M` -- matrix/linear map
- `X` -- convex set

### Notes

This type is parametric in the elements of the linear map, `NM`, which is
independent of the numeric type of the wrapped set (`N`).
Typically `NM = N`, but there may be exceptions, e.g., if `NM` is an interval
that holds numbers of type `N`, where `N` is a floating point number type such
as `Float64`.

### Examples

For the examples we create a ``3×2`` matrix and two unit squares, one of them
being two-dimensional and the other one being one-dimensional.

```jldoctest constructors
julia> A = [1 2; 1 3; 1 4]; X = BallInf([0, 0], 1); Y = BallInf([0], 1);
```

The function ``*`` can be used as an alias to construct a `LinearMap` object.

```jldoctest constructors
julia> lm = LinearMap(A, X)
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([1 2; 1 3; 1 4], BallInf{Int64}([0, 0], 1))

julia> lm2 = A * X
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([1 2; 1 3; 1 4], BallInf{Int64}([0, 0], 1))

julia> lm == lm2
true
```

For convenience, `A` does not need to be a matrix but we also allow to use
vectors (interpreted as an ``n×1`` matrix) and `UniformScaling`s resp. scalars
(interpreted as a scaling, i.e., a scaled identity matrix).
Scaling by ``1`` is ignored.

```jldoctest constructors
julia> using LinearAlgebra: I

julia> [2, 3] * Y
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([2; 3], BallInf{Int64}([0], 1))

julia> lm3 = 2 * X
LinearMap{Int64,BallInf{Int64},Int64,SparseArrays.SparseMatrixCSC{Int64,Int64}}(
  [1, 1]  =  2
  [2, 2]  =  2, BallInf{Int64}([0, 0], 1))

julia> 2I * X == lm3
true

julia> 1I * X == X
true
```

Applying a linear map to a `LinearMap` object combines the two maps into a
single `LinearMap` instance.
Again we can make use of the conversion for convenience.

```jldoctest constructors
julia> B = transpose(A); B * lm
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([3 9; 9 29], BallInf{Int64}([0, 0], 1))

julia> B = [3, 4, 5]; B * lm
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([12 38], BallInf{Int64}([0, 0], 1))

julia> B = 2; B * lm
LinearMap{Int64,BallInf{Int64},Int64,Array{Int64,2}}([2 4; 2 6; 2 8], BallInf{Int64}([0, 0], 1))
```

The application of a `LinearMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest constructors
julia> A * ZeroSet{Int}(2)
ZeroSet{Int64}(3)

julia> A * EmptySet{Int}(2)
EmptySet{Int64}(2)
```
"""
struct LinearMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}} <: LazySet{N}
    M::MAT
    X::S

    # default constructor with dimension match check
    function LinearMap(M::MAT, X::S) where {N<:Real, S<:LazySet{N}, NM,
                                            MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(X))"
        return new{N, S, NM, MAT}(M, X)
    end
end

isoperationtype(::Type{<:LinearMap}) = true
isconvextype(::Type{<:LinearMap{N, S}}) where {N, S} = isconvextype(S)

"""
```
    *(map::Union{AbstractMatrix, UniformScaling, AbstractVector, Real}, X::LazySet)
```

Alias to create a `LinearMap` object.

### Input

- `map` -- linear map
- `X`   -- convex set

### Output

A lazy linear map, i.e., a `LinearMap` instance.
"""
function *(map::Union{AbstractMatrix, UniformScaling, AbstractVector, Real}, X::LazySet)
    return LinearMap(map, X)
end

# convenience constructor from a vector
function LinearMap(v::AbstractVector, X::LazySet)
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
function LinearMap(M::UniformScaling{N}, X::LazySet) where {N<:Real}
    if M.λ == one(N)
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

# disambiguations
function LinearMap(v::AbstractVector, lm::LinearMap)
    return invoke(LinearMap, Tuple{AbstractVector, LazySet}, v, lm)
end

function LinearMap(α::Real, lm::LinearMap)
    return invoke(LinearMap, Tuple{Real, LazySet}, α, lm)
end

# more efficient version
function LinearMap(M::UniformScaling{N}, lm::LinearMap) where {N<:Real}
    if M.λ == one(N)
        return lm
    end
    return LinearMap(M.λ * lm.M, lm.X)
end

# ZeroSet is "almost absorbing" for LinearMap (only the dimension changes)
function LinearMap(M::AbstractMatrix{N}, Z::ZeroSet{N}) where {N<:Real}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

# EmptySet is absorbing for LinearMap
function LinearMap(M::AbstractMatrix{N}, ∅::EmptySet{N}) where {N<:Real}
    return ∅
end

"""
    dim(lm::LinearMap)

Return the dimension of a linear map.

### Input

- `lm` -- linear map

### Output

The ambient dimension of the linear map.
"""
function dim(lm::LinearMap)
    return size(lm.M, 1)
end

"""
    σ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}

Return the support vector of the linear map.

### Input

- `d`  -- direction
- `lm` -- linear map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it follows
that ``σ(d, L) = M⋅σ(M^T d, S)`` for any direction ``d``.
"""
function σ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}
    return lm.M * σ(_At_mul_B(lm.M, d), lm.X)
end

"""
    ρ(d::AbstractVector{N}, lm::LinearMap{N}; kwargs...) where {N<:Real}

Return the support function of the linear map.

### Input

- `d`      -- direction
- `lm`     -- linear map
- `kwargs` -- additional arguments that are passed to the support function
              algorithm

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it follows
that ``ρ(d, L) = ρ(M^T d, S)`` for any direction ``d``.
"""
function ρ(d::AbstractVector{N}, lm::LinearMap{N}; kwargs...) where {N<:Real}
    return ρ(_At_mul_B(lm.M, d), lm.X; kwargs...)
end

"""
    isbounded(lm::LinearMap; cond_tol::Number=DEFAULT_COND_TOL)

Determine whether a linear map is bounded.

### Input

- `lm`       -- linear map
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)

### Output

`true` iff the linear map is bounded.

### Algorithm

We first check if the matrix is zero or the wrapped set is bounded.
If not, we perform a sufficient check whether the matrix is invertible.
If the matrix is invertible, then the map being bounded is equivalent to the
wrapped set being bounded, and hence the map is unbounded.
Otherwise, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(lm::LinearMap; cond_tol::Number=DEFAULT_COND_TOL)
    if iszero(lm.M) || isbounded(lm.X)
        return true
    end
    if isinvertible(lm.M; cond_tol=cond_tol)
        return false
    end
    return isbounded_unit_dimensions(lm)
end

"""
    ∈(x::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}

Check whether a given point is contained in a linear map of a convex set.

### Input

- `x`  -- point/vector
- `lm` -- linear map of a convex set

### Output

`true` iff ``x ∈ lm``.

### Algorithm

Note that ``x ∈ M⋅S`` iff ``M^{-1}⋅x ∈ S``.
This implementation does not explicitly invert the matrix, which is why it also
works for non-square matrices.

### Examples

```jldoctest
julia> lm = LinearMap([2.0 0.0; 0.0 1.0], BallInf([1., 1.], 1.));

julia> [5.0, 1.0] ∈ lm
false
julia> [3.0, 1.0] ∈ lm
true
```

An example with non-square matrix:
```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> [0.5, 0.5] ∈ M*B
true
```
"""
function ∈(x::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}
    return lm.M \ x ∈ lm.X
end

"""
    an_element(lm::LinearMap{N})::Vector{N} where {N<:Real}

Return some element of a linear map.

### Input

- `lm` -- linear map

### Output

An element in the linear map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(lm::LinearMap{N})::Vector{N} where {N<:Real}
    return lm.M * an_element(lm.X)
end

"""
    isempty(lm::LinearMap)

Return if a linear map is empty or not.

### Input

- `lm` -- linear map

### Output

`true` iff the wrapped set is empty.
"""
function isempty(lm::LinearMap)
    return isempty(lm.X)
end

"""
    vertices_list(lm::LinearMap{N}; prune::Bool=true)::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a (polyhedral) linear map.

### Input

- `lm` -- linear map
- `prune` -- (optional, default: `true`) if true removes redundant vertices

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polyhedral.
Then the result is just the linear map applied to the vertices of `X`.
"""
function vertices_list(lm::LinearMap{N}; prune::Bool=true)::Vector{Vector{N}} where {N<:Real}
    # for a zero map, the result is just the list containing the origin
    if iszero(lm.M)
        return [zeros(N, dim(lm))]
    end

    # collect low-dimensional vertices lists
    vlist_X = vertices_list(lm.X)

    # create resulting vertices list
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, lm.M * v)
    end

    return prune ? convex_hull(vlist) : vlist
end

"""
    constraints_list(lm::LinearMap{N}) where {N<:Real}

Return the list of constraints of a (polyhedral) linear map.

### Input

- `lm` -- linear map

### Output

The list of constraints of the linear map.

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

We fall back to a concrete set representation and apply `linear_map`.
"""
function constraints_list(lm::LinearMap{N}) where {N<:Real}
    return constraints_list(linear_map(lm.M, lm.X))
end

"""
    linear_map(M::AbstractMatrix{N}, lm::LinearMap{N}) where {N<:Real}

Return the linear map of a lazy linear map.

### Input

- `M`  -- matrix
- `lm` -- linear map

### Output

The polytope representing the linear map of the lazy linear map of a set.
"""
function linear_map(M::AbstractMatrix{N}, lm::LinearMap{N}) where {N<:Real}
     return linear_map(M * lm.M, lm.X)
end

# lazy projection
function Projection(X::LazySet{N}, vars::AbstractVector{Int}) where {N<:Real}
    M = projection_matrix(N, dim(X), vars)
    return LinearMap(M, X)
end
