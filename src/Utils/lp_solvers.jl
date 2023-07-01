const THREAD_FLOAT_LP_SOLVERs = JuMP.Model[]
const THREAD_EXACT_LP_SOLVERs = JuMP.Model[]

@inline default_lp_solver(::Type{T}) where {T} = default_lp_solver(T, Threads.threadid())
@noinline function default_lp_solver(::Type{T}, tid::Int) where {T}
    THREAD_LP_SOLVERs = thread_specific_lp_solvers(T)

    @assert 0 < tid <= length(THREAD_LP_SOLVERs)
    if @inbounds isassigned(THREAD_LP_SOLVERs, tid)
        @inbounds LP = THREAD_LP_SOLVERs[tid]
    else
        LP = Model(default_lp_solver_factory(T))
        @inbounds THREAD_LP_SOLVERs[tid] = LP
    end
    return LP
end

function init_lp_solvers()
    # Code is heavily inspired by the global state of RNGs in Random.
    # That code contains the following comment, although it is unclear
    # what this pertains to (most likely empty!).
    # "it ensures that we didn't save a bad object"
    resize!(empty!(THREAD_FLOAT_LP_SOLVERs), Threads.nthreads())
    resize!(empty!(THREAD_EXACT_LP_SOLVERs), Threads.nthreads())
    return nothing
end

# default LP solver for floating-point numbers
@inline function default_lp_solver_factory(::Type{<:AbstractFloat})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.SIMPLEX))
end

# default LP solver for rational numbers
@inline function default_lp_solver_factory(::Type{<:Rational})
    return JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.EXACT))
end

@inline thread_specific_lp_solvers(::Type{<:AbstractFloat}) = THREAD_FLOAT_LP_SOLVERs
@inline thread_specific_lp_solvers(::Type{<:Rational}) = THREAD_EXACT_LP_SOLVERs

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
