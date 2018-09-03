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

# constructor for an HPolytope with no constraints
HPolytope{N}() where {N<:Real} = HPolytope{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolytope with no constraints of type Float64
HPolytope() = HPolytope{Float64}()

# conversion constructor
HPolytope(S::LazySet) = convert(HPolytope, S)

# constructor for an HPolytope from a simple H-representation
function HPolytope(A::Matrix{N}, b::Vector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[]
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
    σ(d::AbstractVector{N}, P::HPolytope{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, P::HPolytope{N}) where {N<:Real}
    c = -d
    n = length(constraints_list(P))
    @assert n > 0 "the polytope has no constraints"
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
function addconstraint!(P::HPolytope{N},
                        constraint::LinearConstraint{N})::Nothing where {N<:Real}
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
function tosimplehrep(P::HPolytope{N}) where {N}
    A = hcat([ci.a for ci in P.constraints]...)'
    b = [ci.b for ci in P.constraints]
    return (A, b)
end

function load_polyhedra_hpolytope() # function to be loaded by Requires
return quote

using CDDLib # default backend
import Polyhedra:polyhedron, SimpleHRepresentation, SimpleHRepresentation,
                 HRep, VRep,
                 removehredundancy!, removevredundancy!,
                 hreps, vreps,
                 intersect,
                 convexhull,
                 hcartesianproduct,
                 points

export intersection, convex_hull, cartesian_product, vertices_list, tovrep, tohrep

# HPolytope from an HRep
function HPolytope(P::HRep{N, T}, backend=CDDLib.CDDLibrary()) where {N, T}
    constraints = LinearConstraint{T}[]
    for hi in Polyhedra.allhalfspaces(P)
        push!(constraints, HalfSpace(hi.a, hi.β))
    end
    return HPolytope(constraints)
end

"""
    polyhedron(P::HPolytope{N}, [backend]=CDDLib.CDDLibrary()) where {N}

Return an `HRep` polyhedron from `Polyhedra.jl` given a polytope in H-representation.

### Input

- `P`       -- polytope
- `backend` -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
               computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
               for further information

### Output

An `HRep` polyhedron.
"""
function polyhedron(P::HPolytope{N}, backend=CDDLib.CDDLibrary()) where {N}
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
end

"""
    intersection(P1::HPolytope{N}, P2::HPolytope{N};
                 [backend]=CDDLib.CDDLibrary(),
                 [prunefunc]=removehredundancy!)::HPolytope{N} where {N<:Real}

Compute the intersection of two polytopes in H-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information
- `prunefunc` -- (optional, default: `removehredundancy!`) function to post-process
                  the output of `intersect`

### Output

The `HPolytope` obtained by the intersection of `P1` and `P2`.
"""
function intersection(P1::HPolytope{N}, P2::HPolytope{N};
                      backend=CDDLib.CDDLibrary(),
                      prunefunc=removehredundancy!)::HPolytope{N} where {N<:Real}

    P1 = polyhedron(P1, backend)
    P2 = polyhedron(P2, backend)
    Pint = Polyhedra.intersect(P1, P2)
    prunefunc(Pint)
    return HPolytope(Pint)
end

"""
    convex_hull(P1::HPolytope, P2::HPolytope; [backend]=CDDLib.CDDLibrary())

Compute the convex hull of the set union of two polytopes in H-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `HPolytope` obtained by the concrete convex hull of `P1` and `P2`.
"""
function convex_hull(P1::HPolytope, P2::HPolytope; backend=CDDLib.CDDLibrary())
    Pch = convexhull(polyhedron(P1, backend), polyhedron(P2, backend))
    return HPolytope(Pch)
end

"""
    cartesian_product(P1::HPolytope, P2::HPolytope; [backend]=CDDLib.CDDLibrary())

Compute the Cartesian product of two polytopes in H-representaion.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `HPolytope` obtained by the concrete cartesian product of `P1` and `P2`.
"""
function cartesian_product(P1::HPolytope, P2::HPolytope; backend=CDDLib.CDDLibrary())
    Pcp = hcartesianproduct(polyhedron(P1, backend), polyhedron(P2, backend))
    return HPolytope(Pcp)
end

"""
    vertices_list(P::HPolytope{N};
                  [backend]=CDDLib.CDDLibrary(),
                  [prunefunc]=removevredundancy!)::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polytope in constraint representation.

### Input

- `P`         -- polytope in constraint representation
- `backend`   -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
                 computations backend, see Polyhedra's documentation
                 for further information
- `prunefunc` -- (optional, default: `removevredundancy!`) function to post-process
                 the output of `vreps`

### Output

List of vertices.

### Examples

```jldoctest
julia> using Polyhedra

julia> P = HPolytope([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0], fill(1., 4));

julia> constraints_list(P)
4-element Array{HalfSpace{Float64},1}:
 HalfSpace{Float64}([1.0, 0.0], 1.0)
 HalfSpace{Float64}([0.0, 1.0], 1.0)
 HalfSpace{Float64}([-1.0, 0.0], 1.0)
 HalfSpace{Float64}([0.0, -1.0], 1.0)

julia> vertices_list(P)
4-element Array{Array{Float64,1},1}:
 [1.0, -1.0]
 [1.0, 1.0]
 [-1.0, 1.0]
 [-1.0, -1.0]
```
"""
function vertices_list(P::HPolytope{N};
                       backend=CDDLib.CDDLibrary(),
                       prunefunc=removevredundancy!)::Vector{Vector{N}} where {N<:Real}
    P = polyhedron(P, backend)
    prunefunc(P)
    return collect(points(P))
end

"""
    tovrep(P::HPolytope; backend=CDDLib.CDDLibrary())

Transform a polytope in H-representation to a polytope in V-representation.

### Input

- `P`          -- polytope in constraint representation
- `backend`    -- (optional, default: `CDDLib.CDDLibrary()`) the polyhedral
                  computations backend,
                  see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `VPolytope` which is the vertex representation of the given polytope
in constraint representation.
"""
function tovrep(P::HPolytope; backend=CDDLib.CDDLibrary())
    P = polyhedron(P, backend)
    return VPolytope(P)
end

"""
    tohrep(P::HPolytope)

Return a constraint representation of the given polytope in constraint
representation (no-op).

### Input

- `P` -- polytope in constraint representation

### Output

The same polytope instance.
"""
function tohrep(P::HPolytope)
    return P
end

end # quote
end # function load_polyhedra_hpolytope()
