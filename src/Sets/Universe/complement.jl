function complement(U::Universe{N}) where {N}
    require(@__MODULE__, :LazySets; fun_name="complement")

    return EmptySet{N}(dim(U))
end
