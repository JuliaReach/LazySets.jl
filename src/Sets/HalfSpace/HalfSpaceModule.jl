module HalfSpaceModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, LazySet, AbstractLinearMapAlgorithm
import LinearAlgebra
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: ismultiple, nonzero_indices, samedir
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: isapproxzero, _isapprox, _leq
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: an_element, complement, constraints_list, dim,
                        isbounded, isempty, isoperationtype, isuniversal, rand,
                        distance, ∈, permute, project, ρ, σ, translate
@reexport import ..LazySets: constrained_dimensions, normalize, _is_halfspace,
                             _linear_map_hrep_helper
@reexport import ..Base: convert
@reexport using ..API

export HalfSpace, LinearConstraint,
       halfspace_left, halfspace_right,
       iscomplement

include("HalfSpace.jl")

include("an_element.jl")
include("complement.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("rand.jl")
include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("constrained_dimensions.jl")
include("halfspace_left.jl")
include("halfspace_right.jl")
include("iscomplement.jl")
include("ishalfspace.jl")
include("normalize.jl")

include("convert.jl")

"""
    LinearConstraint

Alias for `HalfSpace`
"""
const LinearConstraint = HalfSpace

"""
    is_tighter_same_dir_2D(c1::HalfSpace,
                           c2::HalfSpace;
                           [strict]::Bool=false)

Check if the first of two two-dimensional constraints with equivalent normal
direction is tighter.

### Input

- `c1`     -- first linear constraint
- `c2`     -- second linear constraint
- `strict` -- (optional; default: `false`) check for strictly tighter
              constraints?

### Output

`true` iff the first constraint is tighter.
"""
function is_tighter_same_dir_2D(c1::HalfSpace,
                                c2::HalfSpace;
                                strict::Bool=false)
    @assert dim(c1) == dim(c2) == 2 "the constraints must be two-dimensional"
    @assert samedir(c1.a, c2.a)[1] "the constraints must have the same " *
                                   "normal direction"

    lt = strict ? (<) : (<=)
    if isapproxzero(c1.a[1])
        @assert isapproxzero(c2.a[1])
        return lt(c1.b, c1.a[2] / c2.a[2] * c2.b)
    end
    return lt(c1.b, c1.a[1] / c2.a[1] * c2.b)
end

# TODO: after #2032, #2041 remove use of this function
_normal_Vector(c::HalfSpace) = HalfSpace(convert(Vector, c.a), c.b)
_normal_Vector(C::Vector{<:HalfSpace}) = [_normal_Vector(c) for c in C]
_normal_Vector(P::LazySet) = _normal_Vector(constraints_list(P))

include("init.jl")

end  # module
