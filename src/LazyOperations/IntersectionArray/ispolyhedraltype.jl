function ispolyhedraltype(::Type{<:IntersectionArray{N,S}}) where {N,S}
    return ispolyhedraltype(S)
end
