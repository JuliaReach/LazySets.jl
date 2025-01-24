import .Polyhedra: polyhedron
export polyhedron, triangulate
using .Polyhedra: HRep, VRep,
                  removehredundancy!, removevredundancy!

function default_polyhedra_backend_1d(N::Type{<:Number}, solver=nothing)
    return Polyhedra.IntervalLibrary{N}()
end

function default_polyhedra_backend_nd(N::Type{<:Number},
                                      solver=default_lp_solver_polyhedra(N))
    return Polyhedra.DefaultLibrary{N}(solver)
end

function default_lp_solver_polyhedra(N::Type{<:AbstractFloat};
                                     presolve::Bool=true)
    if presolve
        return JuMP.optimizer_with_attributes(GLPK.Optimizer,
                                              "presolve" => GLPK_ON)
    else
        return JuMP.optimizer_with_attributes(GLPK.Optimizer)
    end
end

function default_lp_solver_polyhedra(N::Type{<:Rational};
                                     presolve::Bool=false)
    if presolve
        return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT),
                                              "presolve" => GLPK_ON)
    else
        return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT))
    end
end

# solver interface
function _is_polyhedra_backend(backend::Polyhedra.Library)
    return true
end

# in v0.8.0, `Polyhedra` renamed the kwarg `ztol` to `tol`
@static if hasmethod(Polyhedra.detecthlinearity, (Polyhedra.HRepresentation, Any), (:ztol,))
    function _removevredundancy!(X; N)
        return removevredundancy!(X; ztol=_ztol(N))
    end
else
    @assert hasmethod(Polyhedra.detecthlinearity, (Polyhedra.HRepresentation, Any), (:tol,))

    function _removevredundancy!(X; N)
        return removevredundancy!(X; tol=_ztol(N))
    end
end

eval(load_polyhedra_mesh())
eval(load_Polyhedra_polyhedron())

include("init_GeometryBasics.jl")
