function isboundedtype(::Type{<:Rectification{N,S}}) where {N,S}
    return isboundedtype(S)
end
