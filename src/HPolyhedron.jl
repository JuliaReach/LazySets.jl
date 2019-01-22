using MathProgBase, GLPKMathProgInterface

import Base: isempty,
             rand,
             convert

export HPolyhedron,
       dim, σ, ∈,
       addconstraint!,
       constraints_list,
       tosimplehrep,
       tohrep, tovrep,
       convex_hull,
       cartesian_product,
       vertices_list,
       singleton_list,
       isempty,
       linear_map,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       constrained_dimensions

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

# constructor for an HPolyhedron from a simple H-representation
function HPolyhedron(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[]
    @inbounds for i in 1:m
        push!(constraints, LinearConstraint(A[i, :], b[i]))
    end
    return HPolyhedron(constraints)
end

HPolyhedron{N}(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real} = HPolyhedron(A, b)

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
    ρ(d::AbstractVector{N}, P::HPoly{N})::N where {N<:Real}

Evaluate the support function of a polyhedron (in H-representation) in a given
direction.

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation

### Output

The support function of the polyhedron.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result is `Inf`.

### Algorithm

This implementation uses `GLPKSolverLP` as linear programming backend.
"""
function ρ(d::AbstractVector{N}, P::HPoly{N})::N where {N<:Real}
    lp, unbounded = σ_helper(d, P)
    if unbounded
        if P isa HPolytope
            error("the support function in direction $(d) is undefined " *
                  "because the polytope is unbounded")
        end
        return N(Inf)
    end
    return dot(d, lp.sol)
end

"""
    σ(d::AbstractVector{N}, P::HPoly{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, P::HPoly{N}) where {N<:Real}
    lp, unbounded = σ_helper(d, P)
    if unbounded
        if P isa HPolytope
            error("the support vector in direction $(d) is undefined because " *
                  "the polytope is unbounded")
        end
        # construct the solution from the solver's ray result
        ray = (lp == nothing) ? d : lp.attrs[:unboundedray]
        res = Vector{N}(undef, length(ray))
        @inbounds for i in 1:length(ray)
            if ray[i] == zero(N)
                res[i] = zero(N)
            elseif ray[i] > zero(N)
                res[i] = N(Inf)
            else
                res[i] = N(-Inf)
            end
        end
        return res
    else
        return lp.sol
    end
end

@inline function _to_minus_vector(d::SparseVector{N}) where {N}
    c = zeros(N, length(d))
    for (ni, i) in enumerate(d.nzind)
        @inbounds c[i] = -d.nzval[ni]
    end
    return c
end

@inline function _to_minus_vector(d::AbstractVector{N}) where {N}
    return convert(Vector{N}, -d)
end

function σ_helper(d::AbstractVector{N}, P::HPoly{N}) where {N<:Real}

    # let c = -d as a Vector, since GLPK doesn't accept sparse vectors (see #1011)
    c = _to_minus_vector(d)

    (A, b) = tosimplehrep(P)
    if length(b) == 0
        unbounded = true
        lp = nothing
    else
        sense = '<'
        l = -Inf
        u = Inf
        solver = GLPKSolverLP()
        lp = linprog(c, A, sense, b, l, u, solver)
        if lp.status == :Unbounded
            unbounded = true
        elseif lp.status == :Infeasible
            error("the support vector is undefined because the polyhedron is " *
                  "empty")
        else
            unbounded = false
        end
    end
    return (lp, unbounded)
end

"""
    isbounded(P::HPolyhedron)::Bool

Determine whether a polyhedron in constraint representation is bounded.

### Input

- `P` -- polyhedron in constraint representation

### Output

`true` iff the polyhedron is bounded.

### Algorithm

We first check if the polyhedron has more than `max(dim(P), 1)` constraints,
which is a necessary condition for boundedness.
If so, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(P::HPolyhedron)::Bool
    if length(P.constraints) <= max(dim(P), 1)
        return false
    end
    return isbounded_unit_dimensions(P)
end

"""
    ∈(x::AbstractVector{N}, P::HPoly{N})::Bool where {N<:Real}

Check whether a given point is contained in a polyhedron in constraint
representation.

### Input

- `x` -- vector with the coordinates of the point
- `P` -- polyhedron in constraint representation

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

"""
    rand(::Type{HPolyhedron}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::HPolyhedron{N}

Create a polyhedron.

### Input

- `HPolyhedron` -- type for dispatch
- `N`           -- (optional, default: `Float64`) numeric type
- `dim`         -- (optional, default: 2) dimension (is ignored)
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A polyhedron.

### Algorithm

We first create a random polytope and then randomly remove some of the
constraints.
"""
function rand(::Type{HPolyhedron};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::HPolyhedron{N}
    rng = reseed(rng, seed)
    P = rand(HPolytope; N=N, dim=dim, rng=rng)
    constraints_P = constraints_list(P)
    constraints_Q = Vector{LinearConstraint{N}}()
    for i in 1:length(constraints_P)
        if rand(Bool)
            push!(constraints_Q, constraints_P[i])
        end
    end
    return HPolyhedron(constraints_Q)
end

"""
    constrained_dimensions(P::HPolyhedron{N})::Vector{Int} where {N<:Real}

Return the indices in which a polyhedron in constraint representation is
constrained.

### Input

- `P` -- polyhedron in constraint representation

### Output

A vector of ascending indices `i` such that the polyhedron is constrained in
dimension `i`.

### Examples

A 2D polyhedron with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(P::HPolyhedron{N})::Vector{Int} where {N<:Real}
    zero_indices = zeros(Int, dim(P))
    for constraint in P.constraints
        for i in constrained_dimensions(constraint)
            zero_indices[i] = i
        end
    end
    return filter(x -> x != 0, zero_indices)
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
                        constraint::LinearConstraint{N}
                       )::Nothing where {N<:Real}
    push!(P.constraints, constraint)
    return nothing
end

"""
    constraints_list(P::HPoly{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polyhedron in H-representation

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

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` are the
offsets.
"""
function tosimplehrep(P::HPoly{N}) where {N<:Real}
    n = length(constraints_list(P))
    if n == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    A = zeros(N, n, dim(P))
    b = zeros(N, n)
    for (i, Pi) in enumerate(constraints_list(P))
        A[i, :] = Pi.a
        b[i] = Pi.b
    end
    return (A, b)
end

"""
    tohrep(P::HPoly{N}) where {N<:Real}

Return a constraint representation of the given polyhedron in constraint
representation (no-op).

### Input

- `P` -- polyhedron in constraint representation

### Output

The same polyhedron instance.
"""
function tohrep(P::HPoly{N}) where {N<:Real}
    return P
end

"""
    remove_redundant_constraints(P::PT;
                                 backend=GLPKSolverLP())::Union{PT, EmptySet{N}} where {N, PT<:HPoly{N}}

Remove the redundant constraints in a polyhedron in H-representation.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `GLPKSolverLP`) the numeric LP solver backend

### Output

A polyhedron equivalent to `P` but with no redundant constraints, or an empty set
if `P` is detected to be empty (which happens if the constraints are infeasible).

### Algorithm

See [`remove_redundant_constraints!`](@ref) for details.
"""
function remove_redundant_constraints(P::PT;
                                      backend=GLPKSolverLP())::Union{PT, EmptySet{N}} where {N, PT<:HPoly{N}}
    Pred = copy(P)
    if remove_redundant_constraints!(Pred, backend=backend)
        return Pred
    else # the polyhedron P is empty
        return EmptySet{N}()
    end
end

"""
    remove_redundant_constraints!(P::PT;
                                  backend=GLPKSolverLP())::Bool where {N, PT<:HPoly{N}}

Remove the redundant constraints in a polyhedron in H-representation; the polyhedron
is updated in-place.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `GLPKSolverLP`) the numeric LP solver backend

### Output

`true` if the method was successful and the polyhedron `P` is modified by
removing its redundant constraints, and `false` if `P` is detected to be empty
(which happens if the constraints are infeasible).

### Algorithm

If the polyhedron `P` has `m` constraints and its dimension is `n`,
this function checks one by one if each of the `m` constraints is
implied by the remaining ones. To check if the `k`-th constraint
is redundant, an LP is formulated.

For details, see [Fukuda's Polyhedra
FAQ](https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node24.html).
"""
function remove_redundant_constraints!(P::PT;
                                       backend=GLPKSolverLP())::Bool where {N, PT<:HPoly{N}}

    A, b = tosimplehrep(P)
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

    deleteat!(P.constraints, setdiff(1:m, non_redundant_indices))
    return true
end

"""
    linear_map(M::AbstractMatrix{N}, P::PT; [cond_tol=DEFAULT_COND_TOL]::Number)
        where {N<:Real, PT<:HPoly{N}}

Concrete linear map of a polyhedron in constraint representation.

### Input

- `M`        -- matrix
- `P`        -- polyhedron in constraint representation
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)

### Output

A polyhedron of the same type as the input (`PT`).

### Algorithm

If the matrix ``M`` is invertible (which we check with a sufficient condition),
then ``y = M x`` implies ``x = \\text{inv}(M) y`` and we transform the
constraint system ``A x ≤ b`` to ``A \\text{inv}(M) y ≤ b``.
"""
function linear_map(M::AbstractMatrix{N},
                    P::PT;
                    cond_tol::Number=DEFAULT_COND_TOL
                   ) where {N<:Real, PT<:HPoly{N}}
    if !isinvertible(M; cond_tol=cond_tol)
        if P isa HPolyhedron
            error("linear maps for polyhedra need to be invertible")
        end
        # use the implementation for general polytopes
        return invoke(linear_map, Tuple{typeof(M), AbstractPolytope{N}}, M, P)
    end
    # matrix is invertible
    invM = inv(M)
    constraints = Vector{LinearConstraint{N}}(undef, length(constraints_list(P)))
    @inbounds for (i, c) in enumerate(constraints_list(P))
        constraints[i] = LinearConstraint(vec(c.a' * invM), c.b)
    end
    return PT(constraints)
end

# ========================================================
# External methods that require Polyhedra.jl to be loaded
# ========================================================

"""
    convex_hull(P1::HPoly{N}, P2::HPoly{N};
               [backend]=default_polyhedra_backend(P1, N)) where {N}

Compute the convex hull of the set union of two polyhedra in H-representation.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `default_polyhedra_backend(P1, N)`)
                  the polyhedral computations backend

### Output

The `HPolyhedron` (resp. `HPolytope`) obtained by the concrete convex hull of
`P1` and `P2`.

### Notes

For performance reasons, it is suggested to use the `CDDLib.Library()` backend
for the `convex_hull`.

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function convex_hull(P1::HPoly{N},
                     P2::HPoly{N};
                     backend=default_polyhedra_backend(P1, N)) where {N}
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `convex_hull` needs " *
                                               "the package 'Polyhedra' to be loaded"
    Pch = convexhull(polyhedron(P1; backend=backend), polyhedron(P2; backend=backend))
    removehredundancy!(Pch)
    return convert(typeof(P1), Pch)
end

"""
    cartesian_product(P1::HPoly{N}, P2::HPoly{N};
                      [backend]=default_polyhedra_backend(P1, N)
                     ) where {N<:Real}

Compute the Cartesian product of two polyhedra in H-representaion.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `default_polyhedra_backend(P1, N)`)
                  the polyhedral computations backend

### Output

The polyhedron obtained by the concrete cartesian product of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function cartesian_product(P1::HPoly{N},
                           P2::HPoly{N};
                           backend=default_polyhedra_backend(P1, N)
                          ) where {N<:Real}
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `cartesian_product` " *
                                               "needs the package 'Polyhedra' to be loaded"
    Pcp = hcartesianproduct(polyhedron(P1; backend=backend), polyhedron(P2; backend=backend))
    return convert(typeof(P1), Pcp)
end

"""
    tovrep(P::HPoly{N};
          [backend]=default_polyhedra_backend(P, N)) where {N<:Real}

Transform a polyhedron in H-representation to a polytope in V-representation.

### Input

- `P`       -- polyhedron in constraint representation
- `backend` -- (optional, default: `default_polyhedra_backend(P, N)`) the
               backend for polyhedral computations

### Output

The `VPolytope` which is the vertex representation of the given polyhedron
in constraint representation.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1).
"""
function tovrep(P::HPoly{N};
                backend=default_polyhedra_backend(P, N)) where {N<:Real}
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `tovrep` needs " *
                                        "the package 'Polyhedra' to be loaded"
    P = polyhedron(P; backend=backend)
    return VPolytope(P)
end

"""
    vertices_list(P::HPolyhedron{N}) where {N<:Real}

Return the list of vertices of a polyhedron in constraint representation.

### Input

- `P` -- polyhedron in constraint representation

### Output

This function returns an error because the polyhedron is possibly unbounded.
If `P` is known to be bounded, try converting to `HPolytope` first:

```jldoctest
julia> P = HPolyhedron([HalfSpace([1.0, 0.0], 1.0),
                        HalfSpace([0.0, 1.0], 1.0),
                        HalfSpace([-1.0, 0.0], 1.0),
                        HalfSpace([0.0, -1.0], 1.0)]);

julia> P_as_polytope = convert(HPolytope, P);
```
"""
function vertices_list(P::HPolyhedron{N}) where {N<:Real}
    throw(ArgumentError("the list of vertices of a (possibly unbounded) polyhedron is not defined; " *
          "if the polyhedron is bounded, try converting to `HPolytope` first"))
end

"""
    singleton_list(P::HPolyhedron{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polyhedron in H-representation as a list of singletons.

### Input

- `P` -- polytope in constraint representation

### Output

This function returns an error because the polyhedron is possibly unbounded.
If `P` is known to be bounded, try converting to `HPolytope` first.
"""
function singleton_list(P::HPolyhedron{N}) where {N<:Real}
    throw(ArgumentError("the list of singletons of a (possibly unbounded) polyhedron is not defined; " *
          "if the polyhedron is bounded, try converting to `HPolytope` first"))
end

"""
   isempty(P::HPoly{N}; [solver]=GLPKSolverLP())::Bool where {N<:Real}

Determine whether a polyhedron is empty.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `default_polyhedra_backend(P, N)`)
               the polyhedral computations backend
- `solver`  -- (optional, default: `GLPKSolverLP()`) LP solver backend

### Output

`true` if and only if the constraints are inconsistent.

### Algorithm

This function uses `Polyhedra.isempty` which evaluates the feasibility of the
LP whose feasible set is determined by the set of constraints and whose
objective function is zero.

### Notes

This implementation uses `GLPKSolverLP` as linear programming backend by
default.
"""
function isempty(P::HPoly{N};
                 backend=default_polyhedra_backend(P, N),
                 solver=GLPKSolverLP())::Bool where {N<:Real}
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `isempty` needs the " *
                                        "package 'Polyhedra' to be loaded"
    return Polyhedra.isempty(polyhedron(P; backend=backend), solver)
end

convert(::Type{HPolytope}, P::HPolyhedron{N}) where {N<:Real} =
    HPolytope{N}(copy(constraints_list(P)))
convert(::Type{HPolyhedron}, P::HPolytope{N}) where {N<:Real} =
    HPolyhedron{N}(copy(constraints_list(P)))

# ==========================================
# Lower level methods that use Polyhedra.jl
# ==========================================

function load_polyhedra_hpolyhedron() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

@static if VERSION < v"0.7-"
    function convert(::Type{HPolyhedron{N}}, P::HRep{T, N}) where {T, N}
        constraints = LinearConstraint{N}[]
        for hi in Polyhedra.allhalfspaces(P)
            push!(constraints, HalfSpace(hi.a, hi.β))
        end
        return HPolyhedron(constraints)
    end

    function convert(::Type{HPolyhedron}, P::HRep{T, N}) where {T, N}
        return convert(HPolyhedron{N}, P)
    end

    """
        HPolyhedron(P::HRep{T, N}, backend=nothing) where {T, N}

    Return a polyhedron in H-representation given a `HRep` polyhedron
    from `Polyhedra.jl`.

    ### Input

    - `P` -- `HRep` polyhedron

    ### Output

    An `HPolyhedron`.
    """
    function HPolyhedron(P::HRep{T, N}) where {T, N}
        convert(HPolyhedron{N}, P)
    end

else
    function convert(::Type{HPolyhedron{N}}, P::HRep{N}) where {N}
        constraints = LinearConstraint{N}[]
        for hi in Polyhedra.allhalfspaces(P)
            push!(constraints, HalfSpace(hi.a, hi.β))
        end
        return HPolyhedron(constraints)
    end

    function convert(::Type{HPolyhedron}, P::HRep{N}) where {N}
        return convert(HPolyhedron{N}, P)
    end

    function HPolyhedron(P::HRep{N}) where {N}
        convert(HPolyhedron{N}, P)
    end
end # if VERSION < v"0.7-"

"""
    polyhedron(P::HPoly{N};
               [backend]=default_polyhedra_backend(P, N)) where {N<:Real}

Return an `HRep` polyhedron from `Polyhedra.jl` given a polytope in
H-representation.

### Input

- `P`       -- polytope
- `backend` -- (optional, default: call `default_polyhedra_backend(P, N)`)
                the polyhedral computations backend

### Output

An `HRep` polyhedron.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function polyhedron(P::HPoly{N};
                    backend=default_polyhedra_backend(P, N)) where {N<:Real}
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
end

end # quote
end # function load_polyhedra_hpolyhedron()
