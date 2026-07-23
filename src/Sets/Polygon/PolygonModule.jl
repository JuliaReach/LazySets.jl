module PolygonModule

using Reexport: @reexport

using ..LazySets: LazySet

@reexport import ..API: an_element, dim, isconvextype, isbounded, isboundedtype,
                        isempty, isoperationtype, isuniversal, scale, scale!
@reexport using ..API

export Polygon

include("Polygon.jl")

include("an_element.jl")
# include("convex_hull.jl")
include("dim.jl")
include("isconvextype.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
# include("in.jl")
include("scale.jl")
# include("support_function.jl")
# include("support_vector.jl")

end  # module
