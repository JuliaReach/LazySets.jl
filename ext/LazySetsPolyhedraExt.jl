module LazySetsPolyhedraExt

using GLPK: EXACT, Optimizer
using LazySets: GLPK_ON, LazySet, LinearMapElimination, default_cddlib_backend,
                default_polyhedra_backend, isboundedtype, ispolytopic,
                tosimplehrep, _witness_result_empty
using LazySets.HalfSpaceModule: HalfSpace
using LazySets.HPolyhedronModule: HPoly
using LinearAlgebra: I, norm
using JuMP: optimizer_with_attributes
import Polyhedra
using Polyhedra: BlockElimination, DefaultLibrary, HRep, HRepresentation,
                 IntervalLibrary, Library, Polyhedron, allhalfspaces,
                 chebyshevcenter, detecthlinearity, hcartesianproduct, hrep,
                 removehredundancy!, removevredundancy!, volume,
                 hvectortype  # NOTE: this is an internal function
using ReachabilityBase.Comparison: isapproxzero, _ztol
import LazySets: chebyshev_center_radius, default_polyhedra_backend_1d,
                 default_polyhedra_backend_nd, default_lp_solver_polyhedra,
                 polyhedron, _area_Polyhedra, _backend_solver_nd,
                 _cartesian_product_hrep_polyhedra, _get_elimination_instance,
                 _isempty_polyhedron_polyhedra, _is_polyhedra_backend,
                 _minkowski_sum_hrep_preprocess, _removehredundancy!,
                 _removevredundancy!
import LazySets.HPolyhedronModule: _convex_hull

include("Polyhedra/EmptySetExt.jl")
include("Polyhedra/HPolyhedronExt.jl")
include("Polyhedra/HPolytopeExt.jl")
include("Polyhedra/VPolytopeExt.jl")

# in v0.8.0, `Polyhedra` renamed the kwarg `ztol` to `tol`
@static if hasmethod(detecthlinearity, (HRepresentation, Any), (:ztol,))
    function _removevredundancy!(X; N)
        return removevredundancy!(X; ztol=_ztol(N))
    end
else
    @assert hasmethod(detecthlinearity,
                      (HRepresentation, Any),
                      (:tol,)) "there should be a `detecthlinearity` method with `:tol` argument"

    function _removevredundancy!(X; N)
        return removevredundancy!(X; tol=_ztol(N))
    end
end

function default_polyhedra_backend_1d(N::Type{<:Number}, solver=nothing)
    return IntervalLibrary{N}()
end

function default_polyhedra_backend_nd(N::Type{<:Number},
                                      solver=default_lp_solver_polyhedra(N))
    return DefaultLibrary{N}(solver)
end

function default_lp_solver_polyhedra(N::Type{<:AbstractFloat}; presolve::Bool=true)
    if presolve
        return optimizer_with_attributes(Optimizer, "presolve" => GLPK_ON)
    else
        return optimizer_with_attributes(Optimizer)
    end
end

function default_lp_solver_polyhedra(N::Type{<:Rational}; presolve::Bool=false)
    if presolve
        return optimizer_with_attributes(() -> Optimizer(; method=EXACT), "presolve" => GLPK_ON)
    else
        return optimizer_with_attributes(() -> Optimizer(; method=EXACT))
    end
end

# solver interface
function _is_polyhedra_backend(backend::Library)  # NOTE: this is an internal function
    return true
end

"""
    polyhedron(P::LazySet; [backend]=default_polyhedra_backend(P))

Compute a set representation from `Polyhedra.jl`.

### Input

- `P`       -- polyhedral set
- `backend` -- (optional, default: call `default_polyhedra_backend(P)`)
                the polyhedral computations backend

### Output

A set representation in the `Polyhedra` library.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).

### Algorithm

This default implementation uses `tosimplehrep`, which computes the constraint
representation of `P`. Set types preferring the vertex representation should
implement their own method.
"""
function polyhedron(P::LazySet; backend=default_polyhedra_backend(P))
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(hrep(A, b), backend)
end

"""
    chebyshev_center_radius(P::LazySet;
                            [backend]=default_polyhedra_backend(P),
                            [solver]=default_lp_solver_polyhedra(eltype(P); presolve=true))

Compute a [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of a polytopic set.

### Input

- `P`       -- polytopic set
- `backend` -- (optional; default: `default_polyhedra_backend(P)`) the backend
               for polyhedral computations
- `solver`  -- (optional; default:
               `default_lp_solver_polyhedra(N; presolve=true)`) the LP solver
               passed to `Polyhedra`

### Output

The pair `(c, r)` where `c` is a Chebyshev center of `P` and `r` is the radius
of the largest Euclidean ball with center `c` enclosed by `P`.

### Notes

The Chebyshev center is the center of a largest Euclidean ball enclosed by `P`.
In general, the center of such a ball is not unique, but the radius is.

### Algorithm

We call `Polyhedra.chebyshevcenter`.
"""
function chebyshev_center_radius(P::LazySet;
                                 backend=default_polyhedra_backend(P),
                                 solver=default_lp_solver_polyhedra(eltype(P); presolve=true))
    if !ispolytopic(P)
        throw(ArgumentError("can only compute a Chebyshev center for polytopes"))
    end

    Q = polyhedron(P; backend=backend)
    c, r = chebyshevcenter(Q, solver)
    if eltype(P) == Float64  # help with type inference
        c::Vector{Float64}
        r::Float64
    end
    return c, r
end

function _cartesian_product_hrep_polyhedra(P1::PT1, P2::PT2; backend1=nothing,
                                           backend2=nothing) where {PT1<:LazySet,PT2<:LazySet}
    if isnothing(backend1)
        backend1 = default_polyhedra_backend(P1)
    end
    if isnothing(backend2)
        backend2 = default_polyhedra_backend(P2)
    end

    P1′ = polyhedron(P1; backend=backend1)
    P2′ = polyhedron(P2; backend=backend2)
    Pout = hcartesianproduct(P1′, P2′)  # NOTE: this is an internal function

    PT = isboundedtype(PT1) && isboundedtype(PT2) ? HPolytope : HPolyhedron
    return convert(PT, Pout)
end

# common code before calling _minkowski_sum_hrep
function _minkowski_sum_hrep_preprocess(P::LazySet, Q::LazySet, backend, algorithm, prune)
    A, b = tosimplehrep(P)
    C, d = tosimplehrep(Q)
    return _minkowski_sum_hrep(A, b, C, d; backend=backend, algorithm=algorithm, prune=prune)
end

function _backend_solver_nd(N::Type{<:Number})
    backend = default_polyhedra_backend_nd(N)
    solver = default_lp_solver_polyhedra(N)
    return (backend, solver)
end

function _area_Polyhedra(P::LazySet; backend)
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end
    return volume(polyhedron(P; backend=backend))
end

function _get_elimination_instance(N::Type{<:Number}, backend, elimination_method)
    if isnothing(backend)
        backend = default_cddlib_backend(N)
    end
    if isnothing(elimination_method)
        elimination_method = BlockElimination()
    end
    return LinearMapElimination(backend, elimination_method)
end

function _isempty_polyhedron_polyhedra(P::LazySet{N}, witness::Bool=false;
                                       solver=nothing, backend=nothing) where {N}
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    if isnothing(solver)
        result = isempty(polyhedron(P; backend=backend))
    else
        result = isempty(polyhedron(P; backend=backend), solver)
    end

    if result
        return _witness_result_empty(witness, true, N)
    elseif witness
        throw(ArgumentError("witness production is not supported yet"))
    else
        return false
    end
end

function _convex_hull(P1::HPoly, P2::HPoly; backend)
    Pch = convexhull(polyhedron(P1; backend=backend), polyhedron(P2; backend=backend))
    removehredundancy!(Pch)
    return convert(basetype(P1), Pch)
end

# internal helper function for `HPoly`
function _convert_HPoly(T, P::HRep{N}) where {N}
    VT = hvectortype(P)
    constraints = Vector{HalfSpace{N,VT}}()
    for hi in allhalfspaces(P)
        a, b = hi.a, hi.β
        if isapproxzero(norm(a))
            @assert b >= zero(N) "the half-space is inconsistent since it has a zero " *
                                 "normal direction but the constraint is negative"
            continue
        end
        push!(constraints, HalfSpace(a, b))
    end
    return T(constraints)
end

function _removehredundancy!(P::Polyhedron)
    return removehredundancy!(P)
end

end  # module
