# special cases for which the complement is always convex
function isconvextype(::Type{<:Complement{N,<:Union{EmptySet,HalfSpace,Universe}}}) where {N}
    return true
end
