module ZonotopeMDModule

using Reexport, Requires

using ..LazySets: AbstractZonotope
using LinearAlgebra: isdiag, diag
using SparseArrays: blockdiag, sparse, spdiagm

@reexport import ..API: center, isoperationtype, cartesian_product
@reexport import ..LazySets: genmat, ngens

@reexport using ..API

export ZonotopeMD

include("ZonotopeMD.jl")

include("center.jl")
include("isoperationtype.jl")
include("cartesian_product.jl")

include("genmat.jl")
include("ngens.jl")

end  # module
