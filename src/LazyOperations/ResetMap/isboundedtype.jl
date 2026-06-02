function isboundedtype(::Type{<:ResetMap{N,S}}) where {N,S}
    return isboundedtype(S)
end
