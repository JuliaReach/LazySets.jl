module PolygonModule

using Reexport: @reexport

using ..LazySets: LazySet, _plot_recipe_2d_vlist

@reexport import ..API: an_element, dim, isconvextype, isbounded, isboundedtype,
                        isempty, isoperationtype, isuniversal, scale, scale!
import ..LazySets: plot_recipe
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

function plot_recipe(P::Polygon{N}, ε=zero(N)) where {N}
    vlist = P.vertices
    return _plot_recipe_2d_vlist(vlist, N)
end

end  # module
