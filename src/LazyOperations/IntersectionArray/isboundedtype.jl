function isboundedtype(::Type{<:IntersectionArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
