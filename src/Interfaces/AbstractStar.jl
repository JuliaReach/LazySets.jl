"""
    AbstractStar{N} <: LazySet{N}

Abstract supertype for all star set types.
"""
abstract type AbstractStar{N} <: LazySet{N} end

isoperationtype(::Type{<:AbstractStar}) = false
isconvextype(::Type{<:AbstractStar}) = false
