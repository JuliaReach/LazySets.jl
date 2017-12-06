using JuMP, GLPKMathProgInterface

export Polyhedron,
       addconstraint!,
       constraints_list

"""
    Polyhedron{N<:Real} <: LazySet

Type that represents a convex polyhedron in H-representation.

### Fields

- `constraints` -- vector of linear constraints
"""
struct Polyhedron{N<:Real} <: LazySet
    constraints::Vector{LinearConstraint{N}}
end
# constructor for a Polyhedron with no constraints
Polyhedron{N}() where {N<:Real} = Polyhedron{N}(Vector{N}(0))
# constructor for a Polyhedron with no constraints of type Float64
Polyhedron() = Polyhedron{Float64}()

"""
    dim(P::Polyhedron)

Return the dimension of a polyhedron in H-representation.

### Input

- `P`  -- polyhedron in H-representation

### Output

The ambient dimension of the polyhedron in H-representation.
If the polyhedron has no constraints, the result is ``-1``.
"""
function dim(P::Polyhedron)
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    σ(d::AbstractVector{<:Real}, P::Polyhedron)::Vector{<:Real}

Return the support vector of a polyhedron (in H-representation) in a given
direction.

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation

### Output

The support vector in the given direction.

### Algorithm

This implementation uses `GLPKSolverLP` as linear programming backend.
"""
function σ(d::AbstractVector{<:Real}, P::Polyhedron)::Vector{<:Real}
    model = Model(solver=GLPKSolverLP())
    n = length(P.constraints)
    @variable(model, x[1:dim(P)])
    @objective(model, Max, dot(d, x))
    @constraint(model, P[i=1:n],
                dot(P.constraints[i].a, x) <= P.constraints[i].b)
    solve(model)
    return getvalue(x)
end

"""
    addconstraint!(P::Polyhedron{N}, constraint::LinearConstraint{N}) where {N<:Real}

Add a linear constraint to a polyhedron in H-representation.

### Input

- `P`          -- polyhedron in H-representation
- `constraint` -- linear constraint to add

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(P::Polyhedron{N},
                        constraint::LinearConstraint{N}) where {N<:Real}
    push!(P.constraints, constraint)
end

"""
    constraints_list(P::Polyhedron{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polyhedron in H-representation
"""
function constraints_list(P::Polyhedron{N})::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end
