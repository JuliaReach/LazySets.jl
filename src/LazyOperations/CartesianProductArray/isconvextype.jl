function isconvextype(::Type{CartesianProductArray{N,S}}) where {N,S}
    return isconvextype(S)
end
