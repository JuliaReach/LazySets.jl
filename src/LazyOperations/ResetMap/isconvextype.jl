function isconvextype(::Type{ResetMap{N,S}}) where {N,S}
    return isconvextype(S)
end
