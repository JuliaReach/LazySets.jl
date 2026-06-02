function ispolyhedraltype(::Type{<:ConvexHullArray{N,S}}) where {N,S}
    return ispolyhedraltype(S)
end
