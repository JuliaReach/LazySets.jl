function isboundedtype(::Type{<:CartesianProduct{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) && isboundedtype(S2)
end
