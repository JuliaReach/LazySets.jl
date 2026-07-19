using LazySets: LazySet, default_polyhedra_backend, _removevredundancy!
using LazySets.VPolytopeModule: VPolytope
using LazySets.HPolytopeModule: HPolytope
using Polyhedra: VRep, points, removevredundancy, setvrep!, supportssolver,
                 vcartesianproduct, vrep
import Base: convert
import LazySets: polyhedron
import LazySets.VPolytopeModule: _cartesian_product_vrep,
                                 _remove_redundant_vertices, _tohrep

# VPolytope from a VRep
function convert(::Type{VPolytope}, P::VRep)
    vertices = collect(points(P))
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

In this library, `dim` always returns the ambient dimension of the set,
such that a line segment in two dimensions has dimension two. However,
`Polyhedra.dim` will assign a dimension equal to one to a line segment
because it uses a different convention.
"""
function polyhedron(P::VPolytope;
                    backend=default_polyhedra_backend(P),
                    relative_dimension=nothing)
    if isempty(P)
        if isnothing(relative_dimension)
            throw(ArgumentError("the conversion to a `Polyhedra.polyhedron` requires the " *
                                "(relative) dimension of the `VPolytope` to be known, " *
                                "but it cannot be inferred from an empty set; use the " *
                                "keyword argument `relative_dimension`"))
        end
        return Polyhedra.polyhedron(vrep(P.vertices; d=relative_dimension), backend)
    end
    return Polyhedra.polyhedron(vrep(P.vertices), backend)
end

function _cartesian_product_vrep(P1::LazySet, P2::LazySet; backend1=nothing, backend2=nothing)
    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = polyhedron(P1; backend=backend1)
    P2′ = polyhedron(P2; backend=backend2)
    Pout = vcartesianproduct(P1′, P2′)  # NOTE: this is an internal function
    return convert(VPolytope, Pout)
end

function _remove_redundant_vertices(P::VPolytope; backend=nothing, solver=nothing)
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end
    Q = polyhedron(P; backend=backend)
    N = eltype(P)
    if supportssolver(typeof(Q))  # NOTE: this is an internal function
        if isnothing(solver)
            # presolver prints warnings about infeasible solutions (#3226)
            solver = default_lp_solver_polyhedra(N; presolve=false)
        end
        vQ = vrep(Q)
        setvrep!(Q, removevredundancy(vQ, solver))  # NOTE: this is an internal function
    else
        _removevredundancy!(Q; N=N)
    end
    return convert(VPolytope, Q)
end

function _tohrep(P::VPolytope; backend)
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    Q = polyhedron(P; backend=backend)
    R = convert(HPolytope, Q)
    T = _output_type(P)  # help with type inference
    return R::T
end

# internal helper function
function _output_type(::VPolytope)
    return HPolytope{N,Vector{N}} where {N}
end

# internal helper function
function _output_type(::VPolytope{Float64})
    return HPolytope{Float64,Vector{Float64}}
end
