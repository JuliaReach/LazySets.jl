using JuMP.MathOptInterface: AbstractOptimizer, OptimizerWithAttributes
using JuMP: Model, @variable, @objective, @constraint, optimize!,
            termination_status, objective_value, value, primal_status

# solver status
function is_lp_optimal(status)
    return status == JuMP.MathOptInterface.OPTIMAL
end

function is_lp_infeasible(status; strict::Bool=false)
    if status == JuMP.MathOptInterface.INFEASIBLE
        return true
    end
    if strict
        return false
    end
    return status == JuMP.MathOptInterface.INFEASIBLE_OR_UNBOUNDED
end

function is_lp_unbounded(status; strict::Bool=false)
    if status == JuMP.MathOptInterface.DUAL_INFEASIBLE
        return true
    end
    if strict
        return false
    end
    return status == JuMP.MathOptInterface.INFEASIBLE_OR_UNBOUNDED
end

function has_lp_infeasibility_ray(model)
    return primal_status(model) == JuMP.MathOptInterface.INFEASIBILITY_CERTIFICATE
end

# solve a linear program (in the old MathProgBase interface)
function linprog(c, A, sense::Char, b, l::Number, u::Number, solver)
    n = length(c)
    m = length(b)
    return linprog(c, A, fill(sense, m), b, fill(l, n), fill(u, n), solver)
end

function linprog(c, A, sense, b, l, u, solver)
    n = length(c)
    model = Model(solver)
    @variable(model, l[i] <= x[i=1:n] <= u[i])
    @objective(model, Min, c' * x)
    eq_rows, ge_rows, le_rows = sense .== '=', sense .== '>', sense .== '<'
    @constraint(model, A[eq_rows, :] * x .== b[eq_rows])
    @constraint(model, A[ge_rows, :] * x .>= b[ge_rows])
    @constraint(model, A[le_rows, :] * x .<= b[le_rows])
    optimize!(model)
    return (
        status = termination_status(model),
        objval = objective_value(model),
        sol = value.(x),
        model = model
    )
end
