function isconvextype(::Type{MinkowskiSumArray{N,S}}) where {N,S}
    return isconvextype(S)
end
