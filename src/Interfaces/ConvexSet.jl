export ConvexSet

"""
    ConvexSet{N} <: LazySet{N}

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ·x + (1-λ)·y ∈ S``.
"""
abstract type ConvexSet{N} <: LazySet{N} end

isconvextype(X::Type{<:ConvexSet}) = true
