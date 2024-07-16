module HPolyhedronModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, default_lp_solver,
                  default_polyhedra_backend, iscomplement, is_lp_infeasible,
                  is_lp_optimal, is_lp_unbounded, has_lp_infeasibility_ray,
                  linprog, tosimplehrep, _isempty_polyhedron, _normal_Vector
using ..HalfSpaceModule: HalfSpace
using ..HPolytopeModule: HPolytope
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: to_negative_vector
using ReachabilityBase.Basetype: basetype
using ReachabilityBase.Comparison: isapproxzero
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: constraints_list, dim, isempty, isoperationtype, rand,
                        permute, ρ, σ, translate
@reexport import ..LazySets: is_hyperplanar, normalize,
                             remove_redundant_constraints,
                             remove_redundant_constraints!, tohrep, tovrep,
                             addconstraint!
@reexport import ..Base: convert
@reexport using ..API

export HPolyhedron

include("HPolyhedron.jl")

# convenience union type
const HPoly{N} = Union{HPolytope{N},HPolyhedron{N}}

include("constraints_list.jl")
include("dim.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("rand.jl")
include("permute.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("is_hyperplanar.jl")
include("normalize.jl")
include("remove_redundant_constraints.jl")
include("tohrep.jl")
include("tovrep.jl")
include("addconstraint.jl")

function load_polyhedra_hpolyhedron() # function to be loaded by Requires
    return quote
        using .Polyhedra: HRep,
                          polyhedron

        """
             convert(::Type{HPolyhedron}, P::HRep{N}) where {N}

        Convert an `HRep` polyhedron from `Polyhedra.jl` to a polyhedron in constraint
        representation .

        ### Input

        - `HPolyhedron` -- target type
        - `P`           -- `HRep` polyhedron

        ### Output

        An `HPolyhedron`.
        """
        function convert(::Type{HPolyhedron}, P::HRep{N}) where {N}
            VN = Polyhedra.hvectortype(P)
            constraints = Vector{HalfSpace{N,VN}}()
            for hi in Polyhedra.allhalfspaces(P)
                a, b = hi.a, hi.β
                if isapproxzero(norm(a))
                    continue
                end
                push!(constraints, HalfSpace(a, b))
            end
            return HPolyhedron(constraints)
        end

        # convenience conversion method
        function HPolyhedron(P::HRep{N}) where {N}
            return convert(HPolyhedron, P)
        end
    end
end  # quote / load_polyhedra_hpolyhedron()

include("init.jl")

end  # module
