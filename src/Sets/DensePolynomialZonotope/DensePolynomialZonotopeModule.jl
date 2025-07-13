module DensePolynomialZonotopeModule

using Reexport

using ..LazySets: AbstractPolynomialZonotope, _scale_copy_inplace, @validate

@reexport import ..API: center, isoperationtype, linear_map, scale, scale!
@reexport import ..LazySets: ngens_dep, ngens_indep, polynomial_order
@reexport using ..API

export DensePolynomialZonotope

include("DensePolynomialZonotope.jl")

include("center.jl")
include("isoperationtype.jl")
include("linear_map.jl")
include("scale.jl")

include("ngens_dep.jl")
include("ngens_indep.jl")
include("polynomial_order.jl")

end  # module
