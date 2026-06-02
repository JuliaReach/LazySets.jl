function isconvextype(::Type{CachedMinkowskiSumArray{N,S}}) where {N,S}
    return isconvextype(S)
end
