using JuMP, GLPKMathProgInterface

export Polyhedron, addconstraint!, constraints_list

"""
    Polyhedron <: LazySet

Type that represents a convex polyhedron in H-representation.

### Fields

- `constraints` -- a vector of linear constraints
"""
struct Polyhedron <: LazySet
    constraints::Vector{LinearConstraint}
end
Polyhedron() = Polyhedron([])

"""
    dim(P)

Return the ambient dimension of the polyhedron.

### Input

- `P`  -- a polyhedron in H-representation
"""
function dim(P::Polyhedron)::Int64
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    σ(d, P)

Return the support vector of the polyhedron in a given direction.

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation
"""
function σ(d::AbstractVector{Float64}, p::Polyhedron)::Vector{Float64}
    model = Model(solver=GLPKSolverLP())
    n = length(p.constraints)
    @variable(model, x[1:dim(p)])
    @objective(model, Max, dot(d, x))
    @constraint(model, P[i=1:n], dot(p.constraints[i].a, x) <= p.constraints[i].b)
    solve(model)
    return getvalue(x)
end

"""
    addconstraint!(p, c)

Add a linear constraint to a polyhedron.

### Input

- `P`          -- a polyhedron
- `constraint` -- the linear constraint to add

### Notes

It is left to the user to guarantee that the dimension all linear constraints is
the same.
"""
function addconstraint!(P::Polyhedron, c::LinearConstraint)
    push!(P.constraints, c)
end

"""
    constraints_list(P)

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polyhedron in H-representation
"""
function constraints_list(P::Polyhedron)
    return P.constraints
end
