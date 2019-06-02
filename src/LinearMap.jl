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
"""
struct LinearMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}} <: LazySet{N}
    M::MAT
    X::S

    # default constructor with dimension match check
    function LinearMap{N, S, NM, MAT}(M::MAT, X::S) where {N<:Real,
            S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(X))"
        return new{N, S, NM, MAT}(M, X)
    end
end

# convenience constructor without type parameter
LinearMap(M::MAT,
          X::S) where {N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} =
    LinearMap{N, S, NM, MAT}(M, X)

# constructor from a linear map: perform the matrix multiplication immediately
LinearMap(M::MAT, lm::LinearMap{N, S, NM, MAT}) where {N<:Real, S<:LazySet{N},
        NM, MAT<:AbstractMatrix{NM}} =
    LinearMap{N, S, NM, MAT}(M * lm.M, lm.X)

"""
```
    *(M::AbstractMatrix{N}, X::LazySet{N}) where {N<:Real}
```

Return the linear map of a convex set.

### Input

- `M` -- matrix/linear map
- `X` -- convex set

### Output

A lazy linear map, i.e. a `LinearMap` instance.
"""
function *(M::AbstractMatrix{N}, X::LazySet{N}) where {N<:Real}
    return LinearMap(M, X)
end

function *(M::AbstractVector{N}, X::LazySet{N}) where {N<:Real}
    return LinearMap(reshape(M, length(M), 1), X)
end

"""
```
    *(α::N, X::LazySet{N}) where {N<:Real}
```

Return the linear map of a convex set scaled by a given number.

### Input

- `α` -- scalar
- `X` -- convex set

### Output

The lazy linear map of the convex set, `X ↦ M * X`, such that `M` is a sparse
matrix proportional to the identity and whose non-zero entries equal `α`.
"""
function *(α::N, X::LazySet{N}) where {N<:Real}
    n = dim(X)
    return LinearMap(sparse(α * I, n, n), X)
end

"""
```
    *(α::N, lm::LinearMap{N}) where {N<:Real, LM<:LinearMap{N}}
```

Return the linear map scaled by a given value.

### Input

- `α`  -- scalar
- `lm` -- linear map

### Output

The scaled linear map.
"""
function *(α::N, lm::LinearMap{N}) where {N<:Real}
    return LinearMap(α * lm.M, lm.X)
end

"""
```
    *(M::AbstractMatrix{N}, Z::ZeroSet{N})::ZeroSet{N} where {N<:Real}
```

A linear map of a zero set, which is simplified to a zero set (the absorbing
element).

### Input

- `M` -- abstract matrix
- `Z` -- zero set

### Output

The zero set with the output dimension of the linear map.
"""
function *(M::AbstractMatrix{N}, Z::ZeroSet{N})::ZeroSet{N} where {N<:Real}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

"""
    dim(lm::LinearMap)::Int

Return the dimension of a linear map.

### Input

- `lm` -- linear map

### Output

The ambient dimension of the linear map.
"""
function dim(lm::LinearMap)::Int
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
    isbounded(lm::LinearMap; cond_tol::Number=DEFAULT_COND_TOL)::Bool

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
function isbounded(lm::LinearMap; cond_tol::Number=DEFAULT_COND_TOL)::Bool
    if iszero(lm.M) || isbounded(lm.X)
        return true
    end
    if isinvertible(lm.M; cond_tol=cond_tol)
        return false
    end
    return isbounded_unit_dimensions(lm)
end

"""
    ∈(x::AbstractVector{N}, lm::LinearMap{N})::Bool where {N<:Real}

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

julia> ∈([5.0, 1.0], lm)
false
julia> ∈([3.0, 1.0], lm)
true
```

An example with non-square matrix:
```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> ∈([0.5, 0.5], M*B)
true
```
"""
function ∈(x::AbstractVector{N}, lm::LinearMap{N})::Bool where {N<:Real}
    return ∈(lm.M \ x, lm.X)
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
    isempty(lm::LinearMap)::Bool

Return if a linear map is empty or not.

### Input

- `lm` -- linear map

### Output

`true` iff the wrapped set is empty.
"""
function isempty(lm::LinearMap)::Bool
    return isempty(lm.X)
end

"""
    vertices_list(lm::LinearMap{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a (polyhedral) linear map.

### Input

- `lm` -- linear map

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polyhedral.
Then the result is just the linear map applied to the vertices of `X`.
"""
function vertices_list(lm::LinearMap{N})::Vector{Vector{N}} where {N<:Real}
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

    return vlist
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
