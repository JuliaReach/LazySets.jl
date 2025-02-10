module PolygonModule

using Reexport, Requires

using ..LazySets: LazySet, _plot_recipe_2d_vlist
using ReachabilityBase.Comparison: _leq, _geq, _isapprox
using ReachabilityBase.Require: require

@reexport import ..API: convex_hull, dim, isconvextype, isbounded,
                        isboundedtype, isempty, isoperationtype, isuniversal, ∈,
                        ρ, σ
import ..LazySets: plot_recipe
@reexport using ..API

export Polygon

include("Polygon.jl")

include("convex_hull.jl")
include("dim.jl")
include("isconvextype.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("in.jl")
include("support_function.jl")
include("support_vector.jl")

function plot_recipe(P::Polygon{N}, ε=zero(N)) where {N}
    vlist = P.vertices
    return _plot_recipe_2d_vlist(vlist, N)
end

include("init.jl")

end  # module
