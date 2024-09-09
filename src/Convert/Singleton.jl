function convert(::Type{Singleton},
                 cp::CartesianProduct{N,S1,S2}) where {N,S1<:AbstractSingleton,
                                                       S2<:AbstractSingleton}
    return Singleton(vcat(element(first(cp)), element(second(cp))))
end
