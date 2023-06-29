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

# solve a linear program (in the MathOptInterface)
function linprog(c, A, sense, b, l, u, solver_or_model)
    n = length(c)
    model = linprog_model(solver_or_model)

    x = lin_prog_variable!(model, n, l, u)
    lin_prog_constraints!(model, A, sense, b, x)

    @objective(model, Min, c' * x)

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

function lin_prog_variable!(model, n, l::Number, u::Number)
    if isfinite(l) && isfinite(u)
        return @variable(model, x[1:n], lower_bound=l, upper_bound=u)
    elseif isfinite(l)
        return @variable(model, x[1:n], lower_bound=l)
    elseif isfinite(u)
        return @variable(model, x[1:n], upper_bound=u)
    else 
        return @variable(model, x[1:n])
    end
end

function lin_prog_variable!(model, n, l::Vector{T}, u::Vector{T}) where {T}
    return @variable(model, l[i] <= x[i=1:n] <= u[i])
end

function lin_prog_constraints!(model, A, sense::Char, b, x)
    if sense == '='
        @constraint(model, A * x == b)
    elseif sense == '>'
        @constraint(model, A * x >= b)
    elseif sense == '<'
        @constraint(model, A * x <= b)
    else
        throw("Invalid sense ($sense) for constraints in lin_prog.")
    end
end

function lin_prog_constraints!(model, A, sense, b, x)
    eq_rows, ge_rows, le_rows = sense .== '=', sense .== '>', sense .== '<'
    @constraint(model, view(A, eq_rows, :) * x == view(b, eq_rows))
    @constraint(model, view(A, ge_rows, :) * x >= view(b, ge_rows))
    @constraint(model, view(A, le_rows, :) * x <= view(b, le_rows))
end