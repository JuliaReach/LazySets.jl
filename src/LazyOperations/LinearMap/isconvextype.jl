function isconvextype(::Type{<:LinearMap{N,S}}) where {N,S}
    return isconvextype(S)
end
