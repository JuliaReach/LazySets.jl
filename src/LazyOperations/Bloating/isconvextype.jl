function isconvextype(::Type{Bloating{N,S}}) where {N,S}
    return isconvextype(S)
end
