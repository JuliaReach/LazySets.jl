function isboundedtype(::Type{<:Translation{N,S}}) where {N,S}
    return isboundedtype(S)
end
