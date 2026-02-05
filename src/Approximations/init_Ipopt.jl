function default_nln_solver(N::Type{<:Real}=Float64)
    return Ipopt.Optimizer
end

function _underapproximate_box(X::S, solver) where {S<:LazySet}
    if !ispolytopic(X)
        throw(ArgumentError("box underapproximation is only available for convex polytopes"))
    end

    n = dim(X)
    @assert n == 2 "currently only 2D sets are supported"
    P, b = tosimplehrep(X)

    model = Model(solver)
    set_silent(model)
    @variable(model, u[1:n])
    @variable(model, v[1:n])
    @variable(model, x[1:n])
    @variable(model, y[1:n])
    @variable(model, z[1:n])

    # constraints
    @constraint(model, u[2] == 0)
    @constraint(model, v[1] == 0)
    @constraint(model, [i in 1:n], u[i] - y[i] + x[i] == 0)
    @constraint(model, [i in 1:n], v[i] - z[i] + x[i] == 0)
    @constraint(model, P * x .≤ b)
    @constraint(model, P * y .≤ b)
    @constraint(model, P * z .≤ b)
    @constraint(model, P * (y + z - x) .≤ b)

    # maximize area
    @NLobjective(model, Max, u[1] * v[2])
    optimize!(model)

    # convert solution to hyperrectangle
    @inbounds begin
        l = value.(x)
        h = [value(y[1]), value(z[2])]
        if l[1] > h[1]  # sometimes l > h, then just swap them
            tmp = l
            l = h
            h = tmp
        end
    end
    return Hyperrectangle(; low=l, high=h)
end
