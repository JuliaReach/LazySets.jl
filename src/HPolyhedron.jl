using MathProgBase, GLPKMathProgInterface

export HPolyhedron,
       tosimplehrep

"""
    HPolyhedron{N<:Real} <: LazySet{N}

Type that represents a convex polyhedron in H-representation.

### Fields

- `constraints` -- vector of linear constraints
"""
struct HPolyhedron{N<:Real} <: LazySet{N}
    constraints::Vector{LinearConstraint{N}}
end

# constructor for an HPolyhedron with no constraints
HPolyhedron{N}() where {N<:Real} = HPolyhedron{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolyhedron with no constraints of type Float64
HPolyhedron() = HPolyhedron{Float64}()

# constructor for an HPolytope from a simple H-representation
function HPolyhedron(A::Matrix{N}, b::Vector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[]
    @inbounds for i in 1:m
        push!(constraints, LinearConstraint(A[i, :], b[i]))
    end
    return HPolyhedron(constraints)
end

# convenience union type
const HPoly{N} = Union{HPolytope{N}, HPolyhedron{N}}

# --- LazySet interface functions ---

"""
    dim(P::HPoly{N})::Int where {N<:Real}

Return the dimension of a polyhedron in H-representation.

### Input

- `P`  -- polyhedron in H-representation

### Output

The ambient dimension of the polyhedron in H-representation.
If it has no constraints, the result is ``-1``.
"""
function dim(P::HPoly{N})::Int where {N<:Real}
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    σ(d::AbstractVector{N}, P::HPolyhedron{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, P::HPolyhedron{N}) where {N<:Real}
    c = -d
    n = length(constraints_list(P))
    @assert n > 0 "the polyhedron has no constraints"
    A = zeros(N, n, dim(P))
    b = zeros(N, n)
    for (i, Pi) in enumerate(constraints_list(P))
        A[i, :] = Pi.a
        b[i] = Pi.b
    end
    sense = '<'
    l = -Inf
    u = Inf
    solver = GLPKSolverLP()
    lp = linprog(c, A, sense, b, l, u, solver)
    if lp.status == :Unbounded
        error("the support vector in direction $(d) is undefined because " *
              "the polyhedron is unbounded")
    elseif lp.status == :Infeasible
        error("the support vector is undefined because the polyhedron is empty")
    else
        return lp.sol
    end
end

"""
    ∈(x::AbstractVector{N}, P::HPoly{N})::Bool where {N<:Real}

Check whether a given point is contained in a polytope in constraint
representation.

### Input

- `x` -- vector with the coordinates of the point
- `P` -- polytope in constraint representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies on the outside of each hyperplane.
This is equivalent to checking if the point lies in each half-space.
"""
function ∈(x::AbstractVector{N}, P::HPoly{N})::Bool where {N<:Real}
    @assert length(x) == dim(P)

    for c in P.constraints
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end


# ===========================================
# HPolyhedron and HPolytope's shared methods
# ===========================================

"""
    addconstraint!(P::HPoly{N},
                   constraint::LinearConstraint{N})::Nothing where {N<:Real}

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
function addconstraint!(P::HPoly{N},
                        constraint::LinearConstraint{N})::Nothing where {N<:Real}
    push!(P.constraints, constraint)
    return nothing
end

"""
    constraints_list(P::HPoly{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polytope in H-representation

### Output

The list of constraints of the polyhedron.
"""
function constraints_list(P::HPoly{N}
                         )::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end

"""
    tosimplehrep(P::HPoly{N}) where {N}

Return the simple H-representation ``Ax ≤ b`` of a polyhedron.

### Input

- `P` -- polyhedron

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` are the offsets.
"""
function tosimplehrep(P::HPoly{N}) where {N<:Real}
    if length(P.constraints) == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end

    A = hcat([ci.a for ci in P.constraints]...)'
    b = [ci.b for ci in P.constraints]
    return (A, b)
end

