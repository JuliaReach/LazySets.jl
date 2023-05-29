# default semidefinite-programming solver
function default_sdp_solver(::Type{<:Number})
    solver = JuMP.optimizer_with_attributes(() -> SCS.Optimizer())
    JuMP.set_attribute(solver, JuMP.MOI.Silent(), true)
    return solver
end
