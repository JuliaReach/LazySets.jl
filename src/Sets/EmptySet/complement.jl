"""
# Extended help

    complement(∅::EmptySet{N}) where {N}

### Output

The [`Universe`](@ref) of the same dimension.
"""
function complement(∅::EmptySet{N}) where {N}
    require(@__MODULE__, :LazySets; fun_name="complement")

    return Universe{N}(dim(∅))
end
