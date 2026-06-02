function isboundedtype(::Type{<:QuadraticMap{MVT,S}}) where {MVT,S}
    return isboundedtype(S)
end
