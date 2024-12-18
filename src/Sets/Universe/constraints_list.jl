function constraints_list(::Universe{N}) where {N}
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    return HalfSpace{N,Vector{N}}[]
end
