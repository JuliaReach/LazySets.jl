using JuMP: Model, @variable, @objective, @constraint, optimize!,
            termination_status, objective_value, value, primal_status

# solver status
function is_lp_optimal(status)
    return status == JuMP.OPTIMAL
end

function is_lp_infeasible(status; strict::Bool=false)
    if status == JuMP.INFEASIBLE
        return true
    end
    if strict
        return false
    end
    return status == JuMP.INFEASIBLE_OR_UNBOUNDED
end

function is_lp_unbounded(status; strict::Bool=false)
    if status == JuMP.DUAL_INFEASIBLE
        return true
    end
    if strict
        return false
    end
    return status == JuMP.INFEASIBLE_OR_UNBOUNDED
end

function has_lp_infeasibility_ray(model)
    return primal_status(model) == JuMP.INFEASIBILITY_CERTIFICATE
end

# solve a linear program (in the old MathProgBase interface)
function linprog(c, A, sense::Char, b, l::Number, u::Number, solver_or_model)
    n = length(c)
    m = length(b)
    return linprog(c, A, fill(sense, m), b, fill(l, n), fill(u, n), solver_or_model)
end

function linprog(c, A, sense, b, l, u, solver_or_model)
    n = length(c)
    model = linprog_model(solver_or_model)
    @variable(model, l[i] <= x[i=1:n] <= u[i])
    @objective(model, Min, c' * x)
    eq_rows, ge_rows, le_rows = sense .== '=', sense .== '>', sense .== '<'
    @constraint(model, A[eq_rows, :] * x .== b[eq_rows])
    @constraint(model, A[ge_rows, :] * x .>= b[ge_rows])
    @constraint(model, A[le_rows, :] * x .<= b[le_rows])
    optimize!(model)
    return (status = termination_status(model),
            objval = objective_value(model),
            sol    = value.(x),
            model  = model)
end

@inline function linprog_model(model::JuMP.Model)
    Base.empty!(model)
    return model
end
@inline linprog_model(solver) = Model(solver)