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

    function default_lp_solver(N::Type{<:AbstractFloat})
        return JuMP.with_optimizer(GLPK.Optimizer)
    end

    function default_lp_solver(N::Type{<:Rational})
        return JuMP.with_optimizer(GLPK.Optimizer, method=:Exact)
    end
end)

eval(load_polyhedra_hpolytope())
eval(load_polyhedra_hpolyhedron())
eval(load_polyhedra_vpolytope())

eval(initialize_mesh())
