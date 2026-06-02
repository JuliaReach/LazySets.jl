function isconvextype(::Type{Translation{N,S,VN}}) where {N,S,VN}
    return isconvextype(S)
end
