import Base.∈

export constrained_dimensions,
       tosimplehrep,
       remove_redundant_constraints,
       remove_redundant_constraints!

"""
    ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N})::Bool where {N<:Real}

Check whether a given point is contained in a polyhedron.

### Input

- `x` -- point/vector
- `P` -- polyhedron

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies inside each defining half-space.
"""
function ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N})::Bool where {N<:Real}
    @assert length(x) == dim(P) "a $(length(x))-dimensional point cannot be " *
        "an element of a $(dim(P))-dimensional set"

    for c in constraints_list(P)
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end

"""
    constrained_dimensions(P::AbstractPolyhedron)::Vector{Int} where {N<:Real}

Return the indices in which a polyhedron is constrained.

### Input

- `P` -- polyhedron

### Output

A vector of ascending indices `i` such that the polyhedron is constrained in
dimension `i`.

### Examples

A 2D polyhedron with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(P::AbstractPolyhedron)::Vector{Int}
    zero_indices = zeros(Int, dim(P))
    for constraint in constraints_list(P)
        for i in constrained_dimensions(constraint)
            zero_indices[i] = i
        end
    end
    return filter(x -> x != 0, zero_indices)
end

"""
    tosimplehrep(constraints::AbstractVector{LinearConstraint{N}})
        where {N<:Real}

Return the simple H-representation ``Ax ≤ b`` from a list of linear constraints.

### Input

- `constraints` -- a list of linear constraints

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.
"""
function tosimplehrep(constraints::AbstractVector{LinearConstraint{N}}
                     ) where {N<:Real}
    n = length(constraints)
    if n == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    A = zeros(N, n, dim(first(constraints)))
    b = zeros(N, n)
    @inbounds begin
        for (i, Pi) in enumerate(constraints)
            A[i, :] = Pi.a
            b[i] = Pi.b
        end
    end
    return (A, b)
end

"""
     remove_redundant_constraints!(
         constraints::AbstractVector{LinearConstraint{N}};
         [backend]=GLPKSolverLP())::Bool where {N<:Real}

Remove the redundant constraints of a given list of linear constraints; the list
is updated in-place.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `GLPKSolverLP`) numeric LP solver backend

### Output

`true` if the function was successful and the list of constraints `constraints`
is modified by removing the redundant constraints, and `false` only if the
constraints are infeasible.

### Notes

Note that the result may be `true` even if the constraints are infeasible.
For example, ``x ≤ 0 && x ≥ 1`` will return `true` without removing any
constraint.
To check if the constraints are infeasible, use
`isempty(HPolyhedron(constraints)`.

### Algorithm

If there are `m` constraints in `n` dimensions, this function checks one by one
if each of the `m` constraints is implied by the remaining ones.

To check if the `k`-th constraint is redundant, an LP is formulated using the
constraints that have not yet being removed.
If, at an intermediate step, it is detected that a subgroup of the constraints
is infeasible, this function returns `false`.
If the calculation finished successfully, this function returns `true`.

For details, see [Fukuda's Polyhedra
FAQ](https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node24.html).
"""
function remove_redundant_constraints!(constraints::AbstractVector{LinearConstraint{N}};
                                       backend=GLPKSolverLP()
                                      )::Bool where {N<:Real}

    A, b = tosimplehrep(constraints)
    m, n = size(A)
    non_redundant_indices = 1:m

    i = 1 # counter over reduced constraints

    for j in 1:m    # loop over original constraints
        α = A[j, :]
        Ar = A[non_redundant_indices, :]
        br = b[non_redundant_indices]
        br[i] = b[j] + one(N)
        lp = linprog(-α, Ar, '<', br, -Inf, Inf, backend)
        if lp.status == :Infeasible
            # the polyhedron is empty
            return false
        elseif lp.status == :Optimal
            objval = -lp.objval
            if objval <= b[j]
                # the constraint is redundant
                non_redundant_indices = setdiff(non_redundant_indices, j)
            else
                # the constraint is not redundant
                i = i+1
            end
        else
            error("LP is not optimal; the status of the LP is $(lp.status)")
        end
    end

    deleteat!(constraints, setdiff(1:m, non_redundant_indices))
    return true
end

"""
    remove_redundant_constraints(
        constraints::AbstractVector{LinearConstraint{N}};
        backend=GLPKSolverLP())::Union{AbstractVector{LinearConstraint{N}},
                                       EmptySet{N}} where {N<:Real}

Remove the redundant constraints of a given list of linear constraints.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `GLPKSolverLP`) numeric LP solver backend

### Output

The list of constraints with the redundant ones removed, or an empty set if the
constraints are infeasible.

### Algorithm

See
[`remove_redundant_constraints!(::AbstractVector{LinearConstraint{<:Real}})`](@ref)
for details.
"""
function remove_redundant_constraints(constraints::AbstractVector{LinearConstraint{N}};
                                      backend=GLPKSolverLP()
                                     )::Union{AbstractVector{LinearConstraint{N}},
                                                             EmptySet{N}
                                                            } where {N<:Real}
    constraints_copy = copy(constraints)
    if remove_redundant_constraints!(constraints_copy, backend=backend)
        return constraints_copy
    else  # the constraints are infeasible
        return EmptySet{N}()
    end
end
