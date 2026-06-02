function isconvextype(::Type{ExponentialMap{N,S}}) where {N,S}
    return isconvextype(S)
end
