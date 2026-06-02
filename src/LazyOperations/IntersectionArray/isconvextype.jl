function isconvextype(::Type{IntersectionArray{N,S}}) where {N,S}
    return isconvextype(S)
end
