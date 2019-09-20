# default semidefinite-programming solver
function default_SDP_solver(N::Type{<:Real})
    require(:SCS; fun_name="default_SDP_solver")

    return SCS.Optimizer(silent=true)
end
