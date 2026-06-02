function isboundedtype(::Type{<:UnionSetArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
