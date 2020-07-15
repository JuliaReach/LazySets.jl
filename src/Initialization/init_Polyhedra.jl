eval(quote
    import .Polyhedra: polyhedron
    export polyhedron
    using .Polyhedra: HRep, VRep,
                      removehredundancy!, removevredundancy!
    import JuMP, GLPK

    function default_polyhedra_backend(P, N::Type{<:Number})
        if LazySets.dim(P) == 1
            return default_polyhedra_backend_1d(N)
        else
            return default_polyhedra_backend_nd(N)
        end
    end

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
