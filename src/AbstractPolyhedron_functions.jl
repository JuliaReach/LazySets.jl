import Base.∈

export constrained_dimensions,
       tosimplehrep,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       linear_map

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

"""
    linear_map(M::AbstractMatrix{N},
               P::AbstractPolyhedron{N};
               check_invertibility::Bool=true,
               cond_tol::Number=DEFAULT_COND_TOL,
               use_inv::Bool=!issparse(M)
               ) where {N<:Real}

Concrete linear map of a polyhedron in constraint representation.

### Input

- `M` -- matrix
- `P` -- polyhedron in constraint representation
- `check_invertibility` -- (optional, deault: `true`) check if the linear map is
                           invertible, in which case this function uses the matrix
                           inverse; if the invertibility check fails, or if
                           this flag is set to `false`, use the vertex representation
                           to compute the linear map (see below for details)
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)
- `use_inv`  -- (optional, default: `false` if `M` is sparse and `true`
                otherwise) whether to compute the full left division through
                `inv(M)`, or to use the left division for each vector; see below 

### Output

The type of the result is "as close as possible" to the the type of `P`.
To fix the type of the output to something different than the default value,
consider post-processing the result of this function with a call to a suitable
`convert` method.

In particular, the output depends on the type of `P` and the algorithm that was used:

- If the vertex-based approach was used:
    
    - If `P` is a `VPolygon` then the output is a `VPolygon`.
    - If `P` is a `VPolytope` then the output is a `VPolytope`.
    - Otherwise, the output is a `VPolygon` if the dimension of `P` is `2`, or a
      `VPolytope` in other cases.

- If the invertibility criterion was used:

    - The types of `HalfSpace`, `Hyperplane`, `Line` and `AbstractHPolygon` are
      preserved.
    - If `P` is an `AbstractPolytope`, then the output is an `HPolygon` if the
      dimension of `P` is two and an `HPolytope` in other cases.
    - Otherwise, the output is an `HPolyhedron`.

### Algorithm

This function implements two algorithms for the linear map:

- If the matrix ``M`` is invertible (which we check with a sufficient condition),
  then ``y = M x`` implies ``x = \\text{inv}(M) y`` and we transform the
  constraint system ``A x ≤ b`` to ``A \\text{inv}(M) y ≤ b``.
- Otherwise, we transform the polyhedron to vertex representation and apply the map
  to each vertex, returning a polyhedron in vertex representation. 

Note that the vertex representation (second approach) is only available if the
polyhedron is bounded. Hence we check boundedness first.

To switch off the check for invertibility, set the option
`check_invertibility=false`. If `M` is not invertible and the polyhedron
is unbounded, this function returns an exception.

The option `use_inv` lets the user control - in case `M` is invertible - if
the full matrix inverse is computed, or only the left division on the normal
vectors. Note that this helps as a workaround when `M` is a sparse matrix, since
the `inv` function is not available for sparse matrices. In this case, either
use the option `use_inv=false` or convert the type of `M` as in
`linear_map(Matrix(M), P)`.
"""
function linear_map(M::AbstractMatrix{N},
                    P::AbstractPolyhedron{N};
                    check_invertibility::Bool=true,
                    cond_tol::Number=DEFAULT_COND_TOL,
                    use_inv::Bool=!issparse(M)
                   ) where {N<:Real}
    @assert dim(P) == size(M, 2) "a linear map of size $(size(M)) cannot be applied to a set of dimension $(dim(P))"

    if !check_invertibility || !isinvertible(M; cond_tol=cond_tol)
        return _linear_map_vrep(M, P)
    end

    # matrix M is invertible => the normal vectors are vec(c.a' * inv(M))
    return _linear_map_hrep(M, P, use_inv)
end

function _linear_map_vrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    if !isbounded(P)
        throw(ArgumentError("the linear map in vertex representation for an unbounded set is not defined"))
    else
        # since P is bounded, we pass an HPolytope and then convert it to vertex representation
        P = tovrep(HPolytope(constraints_list(P), check_boundedness=false))
    end
    return _linear_map_vrep(M, P)
end

function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                          use_inv::Bool) where {N<:Real}
    constraints = _linear_map_hrep_helper(M, P, use_inv)
    return HPolyhedron(constraints)
end

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                                 use_inv::Bool) where {N<:Real}
    constraints_P = constraints_list(P)
    constraints_MP = similar(constraints_P)
    if use_inv
        invM = inv(M)
        @inbounds for (i, c) in enumerate(constraints_P)
            constraints_MP[i] = LinearConstraint(vec(c.a' * invM), c.b)
        end
    else
        # take left division for each constraint c, transpose(M) \ c.a
        @inbounds for (i, c) in enumerate(constraints_P)
            constraints_MP[i] = LinearConstraint(_At_ldiv_B(M, c.a), c.b)
        end
    end
    return constraints_MP
end
