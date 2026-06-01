function isboundedtype(::Type{<:Bloating{N,S}}) where {N,S}
    return isboundedtype(S)
end
