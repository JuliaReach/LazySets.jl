# special cases for which the complement is always polyhedral
function ispolyhedraltype(::Type{<:Complement{N,<:Union{EmptySet,HalfSpace}}}) where {N}
    return true
end
