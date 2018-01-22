using MathProgBase, GLPKMathProgInterface

export HPolytope,
       addconstraint!,
       constraints_list

"""
    HPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation.

### Fields

- `constraints` -- vector of linear constraints

### Note

This type is more appropriately a *polyhedron*, because no check in the constructor is made
that the constraints determine a bounded set from the finite intersection of half-spaces.
This is a running assumption in this type.
"""
struct HPolytope{N<:Real} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N}}
end
# constructor for a HPolytope with no constraints
HPolytope{N}() where {N<:Real} = HPolytope{N}(Vector{N}(0))
# constructor for a HPolytope with no constraints of type Float64
HPolytope() = HPolytope{Float64}()

"""
    dim(P::HPolytope)::Int

Return the dimension of a polytope in H-representation.

### Input

- `P`  -- polytope in H-representation

### Output

The ambient dimension of the polytope in H-representation.
If it has no constraints, the result is ``-1``.
"""
function dim(P::HPolytope)::Int
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    σ(d::AbstractVector{<:Real}, P::HPolytope)::Vector{<:Real}

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
function σ(d::AbstractVector{<:Real}, P::HPolytope)::Vector{<:Real}
    c = -d
    m = length(constraints_list(P))
    A = zeros(m, dim(P))
    b = zeros(m)
    for (i, Pi) in enumerate(constraints_list(P))
        A[i, :] = Pi.a
        b[i] = Pi.b
    end
    sense = '<'
    l = -Inf
    u = Inf
    solver = GLPKSolverLP()
    lp = linprog(c, A, sense, b, l, u, solver)
    return lp.sol
end

"""
    addconstraint!(P::HPolytope{N},
                   constraint::LinearConstraint{N})::Void where {N<:Real}

Add a linear constraint to a polyhedron in H-representation.

### Input

- `P`          -- polyhedron in H-representation
- `constraint` -- linear constraint to add

### Output

Nothing.

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(P::HPolytope{N},
                        constraint::LinearConstraint{N})::Void where {N<:Real}
    push!(P.constraints, constraint)
    return nothing
end

"""
    constraints_list(P::HPolytope{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polytope in H-representation
"""
function constraints_list(P::HPolytope{N}
                         )::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end

@require Polyhedra begin
    using CDDLib
    import Polyhedra:polyhedron, SimpleVRepresentation, SimpleHRepresentation

    #P1 = polyhedron(SimpleVRepresentation(randn(15, 2)), CDDLib.CDDLibrary())
    #P2 = polyhedron(SimpleVRepresentation(randn(15, 2)), CDDLib.CDDLibrary())
    #Pint = intersect(P1, P2)
    #intersect(P1, P2)

    #function minkowsi_sum(P1, P2)
    #    P1 = polyhedron(SimpleVRepresentation(randn(15, 2)), CDDLibrary())
    #end

    #function convex_hull(P1, P2)
    #
    #end
end
