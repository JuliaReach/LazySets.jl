using MathProgBase, GLPKMathProgInterface

export HPolytope,
       addconstraint!,
       constraints_list,
       tosimplehrep

"""
    HPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation.

### Fields

- `constraints` -- vector of linear constraints

### Note

This type is more appropriately a *polyhedron*, because no check in the
constructor is made that the constraints determine a bounded set from the finite
intersection of half-spaces.
This is a running assumption in this type.
"""
struct HPolytope{N<:Real} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N}}
end
# constructor for a HPolytope with no constraints
HPolytope{N}() where {N<:Real} = HPolytope{N}(Vector{N}(0))
# constructor for a HPolytope with no constraints of type Float64
HPolytope() = HPolytope{Float64}()

# constructor for a HPolytope from a simple H-representation
function HPolytope(A::Matrix{N}, b::Vector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[m]
    @inbounds for i in 1:m
        push!(constraints, LinearConstraint(A[i, :], b[i]))
    end
    return HPolytope(constraints)
end

# --- LazySet interface functions ---


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
    ∈(x::AbstractVector{N}, P::HPolytope{N})::Bool where {N<:Real}

Check whether a given 2D point is contained in a polytope in constraint
representation.

### Input

- `x` -- two-dimensional point/vector
- `P` -- polytope in constraint representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies on the outside of each hyperplane.
This is equivalent to checking if the point lies in each half-space.
"""
function ∈(x::AbstractVector{N}, P::HPolytope{N})::Bool where {N<:Real}
    @assert length(x) == dim(P)

    for c in P.constraints
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end


# --- HPolytope functions ---


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

### Output

The list of constraints of the polyhedron.
"""
function constraints_list(P::HPolytope{N}
                         )::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end

"""
    tosimplehrep(P::HPolytope)

Return the simple H-representation ``Ax ≤ b`` of a polytope.

### Input

- `P` -- polytope

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` are the offsets.
"""
function tosimplehrep(P::HPolytope)
    A = hcat([ci.a for ci in P.constraints]...)'
    b = [ci.b for ci in P.constraints]
    return (A, b)
end

@require Polyhedra begin

using CDDLib # default backend
import Polyhedra:polyhedron, SimpleHRepresentation, SimpleVRepresentation,
                 vreps, hreps,
                 intersect

export intersect

#=
function getsimplehrepresentation(P)
    constraints = LinearConstraint{N}[]
    for hi in hreps(P)
        push!(constraints, LinearConstraint(hi.a, hi.β))
    end
    return HPolytope(constraints)
end
=#

"""
intersect(P1::HPolytope{N}, P2::HPolytope{N},
                       backend=CDDLib.CDDLibrary())::HPolytope where {N<Real}

### Input

- `P1`      -- polytope
- `P2`      -- another polytope
- `backend` -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
               computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
               for further information

### Output

The `HPolytope` obtained by the intersection of `P1` and `P2`.
"""
function intersect(P1::HPolytope{N}, P2::HPolytope{N},
                   backend=CDDLib.CDDLibrary())::HPolytope{N} where {N<:Real}
    P1 = polyhedron(SimpleHRepresentation(tosimplehrep(P1)...), backend)
    P2 = polyhedron(SimpleHRepresentation(tosimplehrep(P2)...), backend)
    Pint = intersect(P1, P2)

    constraints = LinearConstraint{N}[]
    for hi in hreps(Pint)
        push!(constraints, LinearConstraint(hi.a, hi.β))
    end
    return HPolytope(constraints)
end

#function minkowsi_sum(P1, P2)
#    P1 = polyhedron(SimpleVRepresentation(randn(15, 2)), CDDLibrary())
#end

#function convex_hull(P1, P2)
#
#end

end
