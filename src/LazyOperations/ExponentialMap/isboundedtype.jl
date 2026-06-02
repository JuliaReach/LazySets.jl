function isboundedtype(::Type{<:ExponentialMap{N,S}}) where {N,S}
    return isboundedtype(S)
end
