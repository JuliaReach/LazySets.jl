function ispolyhedraltype(::Type{<:Intersection{N,S1,S2}}) where {N,S1,S2}
    return ispolyhedraltype(S1) && ispolyhedraltype(S2)
end
