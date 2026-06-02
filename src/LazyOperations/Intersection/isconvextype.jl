function isconvextype(::Type{Intersection{N,S1,S2}}) where {N,S1,S2}
    return isconvextype(S1) && isconvextype(S2)
end
