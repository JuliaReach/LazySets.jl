module HalfSpaceModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, LazySet, AbstractLinearMapAlgorithm,
                  default_lp_solver, is_lp_infeasible, is_lp_optimal, linprog,
                  _witness_result_empty, @validate, @validate_commutative
import LinearAlgebra
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: ismultiple, nonzero_indices, samedir
using ReachabilityBase.Comparison: isapproxzero, _isapprox, _leq
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: an_element, complement, constraints_list, dim,
                        isbounded, isempty, isoperationtype, isuniversal, rand,
                        distance, ∈, permute, project, ρ, σ, translate,
                        isdisjoint
@reexport import ..LazySets: constrained_dimensions, isfeasible, normalize,
                             remove_redundant_constraints,
                             remove_redundant_constraints!, tosimplehrep
import ..LazySets: _ishalfspace, _linear_map_hrep_helper
import ..Base: convert
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
include("isfeasible.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("rand.jl")
include("remove_redundant_constraints.jl")
include("tosimplehrep.jl")
include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("isdisjoint.jl")

include("halfspace_left.jl")
include("halfspace_right.jl")
include("constrained_dimensions.jl")
include("ishalfspace.jl")
include("normalize.jl")
include("iscomplement.jl")

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
_normal_Vector(P::LazySet) = _normal_Vector(constraints_list(P))
_normal_Vector(c::HalfSpace) = HalfSpace(convert(Vector, c.a), c.b)
_normal_Vector(C::Vector{<:HalfSpace{N}}) where {N} = [_normal_Vector(c) for c in C]

function _normal_Vector(C::Vector{<:HalfSpace})
    N = promote_type([eltype(c) for c in C]...)
    return [HalfSpace(convert(Vector{N}, c.a), N(c.b)) for c in C]
end

include("init.jl")

end  # module
