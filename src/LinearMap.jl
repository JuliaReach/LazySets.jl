export LinearMap

"""
    LinearMap <: LazySet

Type that represents a linear transform of a set. This class is a wrapper
around a linear transformation ``M⋅S`` of a set ``S``, such that it
changes the behaviour of the support vector of the new set.

### Fields

- `M`  -- a linear map, which can a be densem matrix, sparse matrix or a subarray object
- `sf` -- a convex set represented by its support function
"""
struct LinearMap{T<:LazySet, N<:Real} <: LazySet
    M::AbstractMatrix{N}
    sf::T

    # default constructor
    LinearMap{T,N}(M::AbstractMatrix{N}, S::T) where {T<:LazySet,N<:Real} =
        new{T,N}(M, S)
    # in case of constructing a linear map from a linear map, the matrix
    # multiplication is performed here
    LinearMap(M::AbstractMatrix{N}, S::LinearMap{T}) where {T<:LazySet,N<:Real} =
        new{T,N}(M * S.M, S.sf)
end
# type-less convenience constructor
LinearMap(M::AbstractMatrix{N}, sf::T) where {T<:LazySet,N<:Real} = LinearMap{T,N}(M, sf)

import Base.*

"""
    M * X

Linear map of a set.

### Input

- ``M`` -- a matrix, which can be dense or sparse
- ``X`` -- a lazy set

### Output

If the matrix is null, a `ZeroSet` is returned; otherwise a lazy linear map.
"""
function *(M::AbstractMatrix, X::LazySet)
    if findfirst(M) != 0
        return LinearMap(M, X)
    else
        # see also DummySet
        return ZeroSet(dim(X))
    end
end

# linear map of a void set (has to be overridden due to polymorphism reasons)
#=
function *(M::AbstractMatrix, sf::VoidSet)
    if dim(sf) == size(M, 2)
        return VoidSet(size(M, 1))
    else
        throw(DimensionMismatch("a VoidSet of dimension " * string(eval(dim(sf))) *
                " cannot be multiplied by a " * string(eval(size(M))) * " matrix"))
    end
end
=#

"""
    dim(lm)

Ambient dimension of the linear map of a set.

It corresponds to the output dimension of the linear map.

### Input

- `lm` -- a linear map
"""
function dim(lm::LinearMap)::Int64
    return size(lm.M, 1)
end

"""
    σ(d, lm)

Support vector of the linear map of a set.

If `S = MB`, where `M` is sa matrix and `B` is a set, it follows that
`σ(d, S) = Mσ(M^T d, B)` for any direction `d`.

### Input

- `d`  -- a direction
- `lm` -- a linear map
"""
function σ(d::AbstractVector{<:Real}, lm::LinearMap)::Vector{<:Real}
    return lm.M * σ(lm.M.' * d, lm.sf)
end

# multiplication of a set by a scalar value
function *(a::Real, sf::LazySet)
    return LinearMap(sparse(a*I, dim(sf)), sf)
end
