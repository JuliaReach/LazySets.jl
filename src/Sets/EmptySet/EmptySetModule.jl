module EmptySetModule

using Reexport

using ..LazySets: ConvexSet, Universe
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator

@reexport import ..API:  an_element, area, complement, diameter, dim, high, ∈,
                         isbounded, isboundedtype, isconvextype, isempty,
                         isoperationtype, isuniversal, linear_map, low, norm,
                         project, radius, rand, rectify, reflect, scale, scale!,
                         ρ, σ, translate, translate!, vertices, vertices_list,
                         volume
@reexport import ..LazySets: chebyshev_center_radius, plot_recipe

export EmptySet, ∅

include("EmptySet.jl")

include("an_element.jl")
include("area.jl")
include("chebyshev_center_radius.jl")
include("complement.jl")
include("diameter.jl")
include("dim.jl")
include("high.jl")
include("in.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("linear_map.jl")
include("low.jl")
include("norm.jl")
include("project.jl")
include("radius.jl")
include("rand.jl")
include("rectify.jl")
include("reflect.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("vertices_list.jl")
include("vertices.jl")
include("volume.jl")

"""
    plot_recipe(∅::EmptySet{N}, [ε]=zero(N)) where {N}

Convert an empty set to a sequence of points for plotting.
In the special case of an empty set, the sequence is empty.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An empty array.
"""
function plot_recipe(::EmptySet{N}, ε=zero(N)) where {N}
    return []
end

end  # module
