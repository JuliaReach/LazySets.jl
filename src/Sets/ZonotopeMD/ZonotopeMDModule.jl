module ZonotopeMDModule

using Reexport, Requires

using ..LazySets: AbstractZonotope, generators
using LinearAlgebra: isdiag, diag
using SparseArrays: blockdiag, sparse, spdiagm
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: find_unique_nonzero_entry
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: center, isoperationtype, rand, cartesian_product
@reexport import ..LazySets: genmat, ngens
import Base: convert
@reexport using ..API

export ZonotopeMD

include("ZonotopeMD.jl")

include("center.jl")
include("isoperationtype.jl")
include("rand.jl")
include("cartesian_product.jl")

include("genmat.jl")
include("ngens.jl")

include("convert.jl")

end  # module
