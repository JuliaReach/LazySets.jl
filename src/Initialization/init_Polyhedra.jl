eval(quote
    import .Polyhedra: polyhedron
    export polyhedron
    using .Polyhedra: HRep, VRep,
                      removehredundancy!, removevredundancy!
    import JuMP, GLPK

    function default_polyhedra_backend(P, N::Type{<:AbstractFloat})
        return Polyhedra.default_library(LazySets.dim(P), Float64)
    end

    function default_polyhedra_backend(P, N::Type{<:Rational})
        return Polyhedra.default_library(LazySets.dim(P), Rational{Int})
    end

    # NOTE: exists in parallel to `default_lp_solver` because we use different
    # interfaces (see #1493)
    function default_lp_solver_polyhedra(N::Type{<:AbstractFloat})
        return JuMP.with_optimizer(GLPK.Optimizer)
    end

    # NOTE: exists in parallel to `default_lp_solver` because we use different
    # interfaces (see #1493)
    function default_lp_solver_polyhedra(N::Type{<:Rational})
        return JuMP.with_optimizer(GLPK.Optimizer, method=GLPK.EXACT)
    end
end)

eval(load_polyhedra_hpolytope())
eval(load_polyhedra_hpolyhedron())
eval(load_polyhedra_vpolytope())

eval(initialize_mesh())
