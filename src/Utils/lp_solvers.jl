export isfeasible

function default_lp_solver(::Type{T}) where {T}
    key = task_local_lp_solver_key(T)
    solver = get!(() -> JuMP.Model(default_lp_solver_factory(T)), task_local_storage(), key)
    return solver
end

# default LP solver for floating-point numbers
@inline function default_lp_solver_factory(::Type{<:AbstractFloat})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.SIMPLEX))
end

# default LP solver for rational numbers
@inline function default_lp_solver_factory(::Union{Type{<:Rational},Type{Int}})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT))
end

@inline task_local_lp_solver_key(::Type{<:AbstractFloat}) = "LAZYSETS_FLOAT_LP_SOLVER"
@inline task_local_lp_solver_key(::Union{Type{<:Rational},Type{Int}}) = "LAZYSETS_EXACT_LP_SOLVER"

# default LP solver given two possibly different numeric types
@inline function default_lp_solver(M::Type{<:Number}, N::Type{<:Number})
    return default_lp_solver(promote_type(M, N))
end

# check for Polyhedra backend (fallback method)
function _is_polyhedra_backend(backend)
    return false
end

# default LP solver for Polyhedra (fallback method)
# NOTE: exists in parallel to `default_lp_solver` because we use different
# interfaces (see #1493)
function default_lp_solver_polyhedra(N; kwargs...)
    require(@__MODULE__, :Polyhedra; fun_name="default_lp_solver_polyhedra")
    return error("no default solver for numeric type $N")
end

"""
    isfeasible(A::AbstractMatrix, b::AbstractVector, [witness]::Bool=false;
               [solver]=nothing)

Check for feasibility of linear constraints given in matrix-vector form.

### Input

- `A`       -- constraints matrix
- `b`       -- constraints vector
- `witness` -- (optional; default: `false`) flag for witness production
- `solver`  -- (optional; default: `nothing`) LP solver

### Output

If `witness` is `false`, the result is a `Bool`.

If `witness` is `true`, the result is a pair `(res, w)` where `res` is a `Bool`
and `w` is a witness point/vector.

### Algorithm

This implementation solves the corresponding feasibility linear program.
"""
function isfeasible(A::AbstractMatrix, b::AbstractVector, witness::Bool=false;
                    solver=nothing)
    N = promote_type(eltype(A), eltype(b))
    # feasibility LP
    lbounds, ubounds = -Inf, Inf
    sense = '<'
    obj = zeros(N, size(A, 2))
    if isnothing(solver)
        solver = default_lp_solver(N)
    end
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    if is_lp_optimal(lp.status)
        return witness ? (true, lp.sol) : true
    elseif is_lp_infeasible(lp.status)
        return _witness_result_empty(witness, false, N)
    end
    return error("LP returned status $(lp.status) unexpectedly")
end
