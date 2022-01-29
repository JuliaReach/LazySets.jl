using JuMP.MathOptInterface: AbstractOptimizer, OptimizerWithAttributes
using JuMP: Model, @variable, @objective, @constraint, optimize!, termination_status, objective_value, value 

# solver statuses
const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.INFEASIBLE_OR_UNBOUNDED

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
        sol = value.(x)
    )
end
