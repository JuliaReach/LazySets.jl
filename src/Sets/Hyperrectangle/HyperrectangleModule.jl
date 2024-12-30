module HyperrectangleModule

using Reexport, Requires

using ..LazySets: AbstractHyperrectangle, _scale_copy_inplace,
                  _ρ_sev_hyperrectangle, _σ_sev_hyperrectangle
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityBase.Comparison: isapproxzero
using ReachabilityBase.Distribution: reseed!
using SparseArrays: SparseVector, findnz, sparse

@reexport import ..API: center, isoperationtype, rand, permute, scale, scale!,
                        ρ, σ, translate, translate!
@reexport import ..LazySets: genmat, radius_hyperrectangle, □
import ..LazySets: _genmat_static
import Base: convert
@reexport using ..API

export Hyperrectangle

include("Hyperrectangle.jl")

include("center.jl")
include("isoperationtype.jl")
include("rand.jl")
include("permute.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("genmat.jl")
include("radius_hyperrectangle.jl")

include("convert.jl")

function □(c::VNC, r::VNR) where {N,VNC<:AbstractVector{N},VNR<:AbstractVector{N}}
    return Hyperrectangle(c, r)
end

include("init.jl")

end  # module
