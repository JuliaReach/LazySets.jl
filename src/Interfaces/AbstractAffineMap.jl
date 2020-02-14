import Base: isempty, ∈

export AbstractAffineMap,
       get_A, get_b, get_X

"""
    AbstractAffineMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Abstract type for affine maps.

### Notes

See [`AffineMap`](@ref) for a standard implementation of this interface.

Every concrete `AbstractAffineMap` must define the following functions:
- `get_A(::AbstractAffineMap)` -- return the linear map
- `get_b(::AbstractAffineMap)` -- return the affine translation vector
- `get_X(::AbstractAffineMap)` -- return the set that the map is applied to
"""
abstract type AbstractAffineMap{N<:Real, S<:LazySet{N}} <: LazySet{N} end

isoperationtype(::Type{<:AbstractAffineMap}) = true
isconvextype(::Type{<:AbstractAffineMap{N, S}}) where {N, S} = isconvextype(S)


# --- AbstractAffineMap interface functions ---


