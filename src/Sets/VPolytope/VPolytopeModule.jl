module VPolytopeModule

using Reexport, Requires

using ..LazySets: AbstractPolytope, LazySet, LinearMapVRep, default_lp_solver,
                  default_lp_solver_polyhedra, default_polyhedra_backend,
                  is_lp_infeasible, is_lp_optimal, linprog
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: projection_matrix
using ReachabilityBase.Comparison: _ztol
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: constraints_list, dim, isoperationtype, rand, reflect,
                        vertices_list, ∈, linear_map, permute, project, scale!,
                        ρ, σ, translate, translate!
@reexport import ..LazySets: remove_redundant_vertices, tohrep, tovrep,
                             _linear_map_vrep
import Base: convert
@reexport using ..API

export VPolytope

include("VPolytope.jl")

include("constraints_list.jl")
include("dim.jl")
include("rand.jl")
include("reflect.jl")
include("vertices_list.jl")
include("in.jl")
include("isoperationtype.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("remove_redundant_vertices.jl")
include("tohrep.jl")
include("tovrep.jl")

include("convert.jl")

function load_polyhedra_vpolytope() # function to be loaded by Requires
    return quote
        import .Polyhedra: polyhedron
        using .Polyhedra: VRep

        # VPolytope from a VRep
        function VPolytope(P::VRep{N}) where {N}
            vertices = collect(Polyhedra.points(P))
            return VPolytope(vertices)
        end

        """
            polyhedron(P::VPolytope;
                       [backend]=default_polyhedra_backend(P),
                       [relative_dimension]=nothing)

        Return a `VRep` polyhedron from `Polyhedra.jl` given a polytope in vertex
        representation.

        ### Input

        - `P`       -- polytope in vertex representation
        - `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
                       backend for polyhedral computations; see [Polyhedra's
                       documentation](https://juliapolyhedra.github.io/) for further
                       information
        - `relative_dimension` -- (default, optional: `nothing`) an integer representing
                                  the (relative) dimension of the polytope; this
                                  argument is mandatory if the polytope is empty

        ### Output

        A `VRep` polyhedron.

        ### Notes

        The *relative dimension* (or just *dimension*) refers to the dimension of the
        set relative to itself, independently of the ambient dimension. For example, a
        point has (relative) dimension zero, and a line segment has (relative) dimension
        one.

        In this library, `LazySets.dim` always returns the ambient dimension of the set,
        such that a line segment in two dimensions has dimension two. However,
        `Polyhedra.dim` will assign a dimension equal to one to a line segment
        because it uses a different convention.
        """
        function polyhedron(P::VPolytope;
                            backend=default_polyhedra_backend(P),
                            relative_dimension=nothing)
            if isempty(P)
                if isnothing(relative_dimension)
                    error("the conversion to a `Polyhedra.polyhedron` requires the " *
                          "(relative) dimension of the `VPolytope` to be known, but it " *
                          "cannot be inferred from an empty set; use the keyword " *
                          "argument `relative_dimension`")
                end
                return polyhedron(Polyhedra.vrep(P.vertices; d=relative_dimension),
                                  backend)
            end
            return polyhedron(Polyhedra.vrep(P.vertices), backend)
        end
    end
end  # quote / load_polyhedra_vpolytope()

include("init.jl")

end  # module
