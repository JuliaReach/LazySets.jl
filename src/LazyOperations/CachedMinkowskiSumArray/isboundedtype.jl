function isboundedtype(::Type{<:CachedMinkowskiSumArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
