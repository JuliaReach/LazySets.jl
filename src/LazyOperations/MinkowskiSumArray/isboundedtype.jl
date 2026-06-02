function isboundedtype(::Type{<:MinkowskiSumArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
