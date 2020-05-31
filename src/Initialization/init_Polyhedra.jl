eval(quote
    import .Polyhedra: polyhedron
    export polyhedron
    using .Polyhedra: HRep, VRep,
                      removehredundancy!, removevredundancy!
    import JuMP, GLPK

    function default_polyhedra_backend(P, N::Type{<:AbstractFloat})
        return Polyhedra.default_library(LazySets.dim(P), N)
    end

    function default_polyhedra_backend(P, N::Type{<:Rational})
        return Polyhedra.default_library(LazySets.dim(P), N)
    end

    function default_lp_solver_polyhedra(N::Type{<:AbstractFloat};
                                         presolve::Bool=false)
        if presolve
            return JuMP.optimizer_with_attributes(GLPK.Optimizer,
                                                  "presolve" => GLPK.ON)
        else
            return JuMP.optimizer_with_attributes(GLPK.Optimizer)
        end
    end

    function default_lp_solver_polyhedra(N::Type{<:Rational};
                                         presolve::Bool=false)
        if presolve
            return JuMP.optimizer_with_attributes(
                () -> GLPK.Optimizer(method=GLPK.EXACT),
                "presolve" => GLPK.ON)
        else
            return JuMP.optimizer_with_attributes(
                () -> GLPK.Optimizer(method=GLPK.EXACT))
        end
    end
end)

eval(load_polyhedra_hpolytope())
eval(load_polyhedra_hpolyhedron())
eval(load_polyhedra_vpolytope())
eval(load_polyhedra_mesh())
