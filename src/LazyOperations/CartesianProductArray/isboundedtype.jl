function isboundedtype(::Type{<:CartesianProductArray{N,S}}) where {N,S}
    return isboundedtype(S)
end
