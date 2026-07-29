import LazySets

using LazySets.EmptySetModule: EmptySet
using LazySets.UniverseModule: Universe
using LazySets: dim
import LazySets.API: complement

"""
# Extended help

    complement(∅::EmptySet)

### Output

The [`Universe`](@ref) of the same dimension.
"""
function complement(∅::EmptySet)
    N = eltype(∅)
    return Universe{N}(dim(∅))
end
