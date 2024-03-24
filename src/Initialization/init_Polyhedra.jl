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

function default_lp_solver_polyhedra(::Type{<:AbstractFloat};
                                     presolve::Bool=false)
    if presolve
        return JuMP.optimizer_with_attributes(GLPK.Optimizer,
                                              "presolve" => GLPK_ON)
    else
        return JuMP.optimizer_with_attributes(GLPK.Optimizer)
    end
end

function default_lp_solver_polyhedra(::Type{<:Rational};
                                     presolve::Bool=false)
    if presolve
        return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT),
                                              "presolve" => GLPK_ON)
    else
        return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT))
    end
end

# solver interface
function _is_polyhedra_backend(::Polyhedra.Library)
    return true
end

eval(load_polyhedra_hpolytope())
eval(load_polyhedra_hpolyhedron())
eval(load_polyhedra_vpolytope())
eval(load_polyhedra_universe())
eval(load_polyhedra_mesh())
eval(load_polyhedra_lazyset())
