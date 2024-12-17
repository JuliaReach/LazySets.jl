"""
# Extended help

    complement(∅::EmptySet)

### Output

The [`Universe`](@ref) of the same dimension.
"""
function complement(∅::EmptySet)
    require(@__MODULE__, :LazySets; fun_name="complement")

    N = eltype(∅)
    return Universe{N}(dim(∅))
end
