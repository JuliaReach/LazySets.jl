function ispolyhedraltype(::Type{<:CartesianProductArray{N,S}}) where {N,S}
    return ispolyhedraltype(S)
end
