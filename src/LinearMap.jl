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
type LinearMap <: LazySet
    M::Union{Matrix{Float64}, SparseMatrixCSC{Float64,Int64}, SubArray}
    sf::LazySet

    # in case of constructing a linear map from a linear map, the matrix
    # multiplication is performed here 
    function LinearMap(M, S)
        if isa(S, LinearMap)
            return new(M * S.M, S.sf)
        else
            return new(M, S)
        end
    end
end

import Base.*

# linear map of a set
function *(M::Union{Matrix, SparseMatrixCSC, SubArray}, sf::LazySet)
    if findfirst(M) != 0
        return LinearMap(M, sf)
    else
        return VoidSet(dim(sf))
    end
end

# linear map of a void set (has to be overridden due to polymorphism reasons)
function *(M::Union{Matrix, SparseMatrixCSC}, sf::VoidSet)
    if dim(sf) == size(M, 2)
        return VoidSet(size(M, 1))
    else
        throw(DimensionMismatch("a VoidSet of dimension " * string(eval(dim(sf))) *
                " cannot be multiplied by a " * string(eval(size(M))) * " matrix"))
    end
end

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
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, lm::LinearMap)::Vector{Float64}
    return lm.M * σ(lm.M.' * d, lm.sf)
end

# multiplication of a set by a scalar value
function *(a::Float64, sf::LazySet)
    return LinearMap(sparse(a*I, dim(sf)), sf)
end

