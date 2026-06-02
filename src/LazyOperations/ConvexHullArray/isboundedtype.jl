function isboundedtype(::Type{<:ConvexHullArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
