module UniverseModule

using Reexport, Requires

using ..LazySets: LazySet, AbstractPolyhedron, default_polyhedra_backend
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator
using ReachabilityBase.Require: require

@reexport import ..API: an_element, complement, constraints, constraints_list,
                        diameter, dim, isbounded, isboundedtype, isempty,
                        isoperationtype, isuniversal, norm, radius, rand,
                        reflect, volume, ∈, permute, project, scale, scale!, ρ,
                        σ, translate, translate!, cartesian_product,
                        convex_hull, intersection
@reexport import ..LazySets: constrained_dimensions, linear_map_inverse,
                             rationalize, tosimplehrep
import Base: copy
@reexport using ..API

export Universe

include("Universe.jl")

include("an_element.jl")
include("complement.jl")
include("constraints.jl")
include("constraints_list.jl")
include("copy.jl")
include("diameter.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("norm.jl")
include("radius.jl")
include("rand.jl")
include("rationalize.jl")
include("reflect.jl")
include("volume.jl")
include("in.jl")
include("permute.jl")
include("project.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("cartesian_product.jl")
include("convex_hull.jl")
include("intersection.jl")

include("constrained_dimensions.jl")
include("linear_map_inverse.jl")
include("tosimplehrep.jl")

include("init.jl")

end  # module
