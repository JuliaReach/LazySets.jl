function isboundedtype(::Type{<:LinearMap{N,S}}) where {N,S}
    return isboundedtype(S)
end
