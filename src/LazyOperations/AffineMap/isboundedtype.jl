function isboundedtype(::Type{<:AffineMap{N,S}}) where {N,S}
    return isboundedtype(S)
end
