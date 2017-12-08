import Base.*

export LinearMap

"""
    LinearMap{S<:LazySet, N<:Real} <: LazySet

Type that represents a linear transformation ``M⋅S`` of a convex set ``S``.

### Fields

- `M`  -- matrix/linear map
- `sf` -- convex set
"""
struct LinearMap{S<:LazySet, N<:Real} <: LazySet
    M::AbstractMatrix{N}
    sf::S
end
# constructor from a linear map: perform the matrix multiplication immediately
LinearMap(M::AbstractMatrix{N}, map::LinearMap{S}) where {S<:LazySet, N<:Real} =
    LinearMap{S,N}(M * map.M, map.sf)

"""
```
    *(M::AbstractMatrix{<:Real}, S::LazySet)
```

Return the linear map of a convex set.

### Input

- `M` -- matrix/linear map
- `S` -- convex set

### Output

The linear map of the convex set.
"""
function *(M::AbstractMatrix{<:Real}, S::LazySet)
    if findfirst(M) != 0
        return LinearMap(M, S)
    else
        return VoidSet(dim(S))
    end
end

"""
```
    *(a::Real, S::LazySet)::LinearMap
```

Return a linear map of a convex set by a scalar value.

### Input

- `a` -- real scalar
- `S` -- convex set

### Output

The linear map of the convex set.
"""
function *(a::Real, S::LazySet)::LinearMap
    return LinearMap(sparse(a*I, dim(S)), S)
end

# linear map of a void set (has to be overridden due to polymorphism reasons)
function *(M::AbstractMatrix, sf::VoidSet)
    if dim(sf) == size(M, 2)
        return VoidSet(size(M, 1))
    else
        throw(DimensionMismatch("a VoidSet of dimension " * string(eval(dim(sf))) *
                " cannot be multiplied by a " * string(eval(size(M))) * " matrix"))
    end
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
    σ(d::AbstractVector{<:Real}, lm::LinearMap)::AbstractVector{<:Real}

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
function σ(d::AbstractVector{<:Real}, lm::LinearMap)::AbstractVector{<:Real}
    return lm.M * σ(lm.M.' * d, lm.sf)
end
