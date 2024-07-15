module HPolytopeModule

using Reexport, Requires

using ..LazySets: AbstractPolytope, HalfSpace, HPolygon,
                  AbstractLinearMapAlgorithm, default_polyhedra_backend,
                  vertices_list_1d, _linear_map_hrep, _normal_Vector
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero, _ztol
using ReachabilityBase.Require: require

@reexport import ..API: isbounded, isoperationtype, rand, vertices_list
@reexport import ..LazySets: _linear_map_hrep_helper, _vertices_list
@reexport import Base: convert
@reexport using ..API

export HPolytope

include("HPolytope.jl")

include("isbounded.jl")
include("isoperationtype.jl")
include("rand.jl")
include("vertices_list.jl")
include("linear_map.jl")

function load_polyhedra_hpolytope() # function to be loaded by Requires
    return quote
        using .Polyhedra: HRep

        function convert(::Type{HPolytope}, P::HRep{N}) where {N}
            VT = Polyhedra.hvectortype(P)
            constraints = Vector{HalfSpace{N,VT}}()
            for hi in Polyhedra.allhalfspaces(P)
                a, b = hi.a, hi.β
                if isapproxzero(norm(a))
                    @assert b >= zero(N) "the half-space is inconsistent since it " *
                                         "has a zero normal direction but the constraint is negative"
                    continue
                end
                push!(constraints, HalfSpace(hi.a, hi.β))
            end
            return HPolytope(constraints)
        end

        """
            HPolytope(P::HRep)

        Return a polytope in constraint representation given an `HRep` polyhedron from
        `Polyhedra.jl`.

        ### Input

        - `P` -- `HRep` polyhedron

        ### Output

        An `HPolytope`.
        """
        function HPolytope(P::HRep)
            return convert(HPolytope, P)
        end
    end
end  # quote / load_polyhedra_hpolytope()

include("init.jl")

end  # module
