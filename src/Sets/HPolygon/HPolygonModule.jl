module HPolygonModule

using Reexport

using ..LazySets: AbstractHPolygon, BINARY_SEARCH_THRESHOLD, addconstraint!,
                  binary_search_constraints, constraints_list, element,
                  isbounded, ⪯, _intersection_line2d
using ..HalfSpaceModule: HalfSpace, _normal_Vector

@reexport import ..API: isoperationtype, σ, translate
@reexport using ..API

export HPolygon

include("HPolygon.jl")

include("isoperationtype.jl")
include("support_vector.jl")
include("translate.jl")

end  # module
