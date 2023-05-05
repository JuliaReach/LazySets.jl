# default LP solver for floating-point numbers
function default_lp_solver(::Type{<:AbstractFloat})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.SIMPLEX))
end

# default LP solver for rational numbers
function default_lp_solver(::Type{<:Rational})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT))
end

# default LP solver given two possibly different numeric types
function default_lp_solver(M::Type{<:Number}, N::Type{<:Number})
    return default_lp_solver(promote_type(M, N))
end

# check for Polyhedra backend (fallback method)
function _is_polyhedra_backend(backend)
    return false
end

# default LP solver for Polyhedra (fallback method)
# NOTE: exists in parallel to `default_lp_solver` because we use different
# interfaces (see #1493)
function default_lp_solver_polyhedra(N, kwargs...)
    require(@__MODULE__, :Polyhedra; fun_name="default_lp_solver_polyhedra")
    return error("no default solver for numeric type $N")
end
