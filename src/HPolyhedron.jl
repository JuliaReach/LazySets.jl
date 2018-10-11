using MathProgBase, GLPKMathProgInterface

export HPolyhedron,
       dim, σ, ∈,
       addconstraint!,
       constraints_list,
       tosimplehrep,
       tohrep, tovrep,
       convex_hull,
       cartesian_product,
       vertices_list
       
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
    c = -d
    (A, b) = tosimplehrep(P)
    if length(b) == 0
        unbounded = true
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
        end
        unbounded = false
    end
    if unbounded
        error("the support vector in direction $(d) is undefined because " *
              "the polyhedron is unbounded")
    else
        return lp.sol
    end
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

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` are the offsets.
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
    tohrep(P::HPoly{N}) where {N}

Return a constraint representation of the given polyhedron in constraint
representation (no-op).

### Input

- `P` -- polyhedron in constraint representation

### Output

The same polyhedron instance.
"""
function tohrep(P::HPoly{N}) where {N}
    return P
end

# ========================================================
# External methods that require Polyhedra.jl to be loaded
# ========================================================

"""
    convex_hull(P1::HPoly{N}, P2::HPoly{N};
               [backend]=default_polyhedra_backend(N)) where {N}

Compute the convex hull of the set union of two polyhedra in H-representation.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`)
                  the polyhedral computations backend

### Output

The `HPolyhedron` (resp. `HPolytope`) obtained by the concrete convex hull of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function convex_hull(P1::HPoly{N}, P2::HPoly{N}; backend=default_polyhedra_backend(N)) where {N}
    @assert isdefined(Main, :Polyhedra) "the function `convex_hull` needs " *
                                        "the package 'Polyhedra' to be loaded"
    Pch = convexhull(polyhedron(P1, backend), polyhedron(P2, backend))
    return convert(typeof(P1), Pch)
end

"""
    cartesian_product(P1::HPOLY, P2::HPOLY;
                      [backend]=default_polyhedra_backend(N)) where {N, HPOLY<:HPolytope{N}}

Compute the Cartesian product of two polyhedra in H-representaion.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`)
                  the polyhedral computations backend

### Output

The polyhedron obtained by the concrete cartesian product of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function cartesian_product(P1::HPoly{N}, P2::HPoly{N};
                          backend=default_polyhedra_backend(N)) where {N}
    @assert isdefined(Main, :Polyhedra) "the function `cartesian_product` needs " *
                                        "the package 'Polyhedra' to be loaded"
    Pcp = hcartesianproduct(polyhedron(P1, backend), polyhedron(P2, backend))
    return convert(typeof(P1), Pcp)
end

"""
    tovrep(P::HPoly{N};
          [backend]=default_polyhedra_backend(N)) where {N}

Transform a polyhedron in H-representation to a polytope in V-representation.

### Input

- `P`          -- polyhedron in constraint representation
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`)
                  the polyhedral computations backend

### Output

The `VPolytope` which is the vertex representation of the given polyhedron
in constraint representation.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function tovrep(P::HPoly{N}; backend=default_polyhedra_backend(N)) where {N}
    @assert isdefined(Main, :Polyhedra) "the function `tovrep` needs " *
                                        "the package 'Polyhedra' to be loaded"
    P = polyhedron(P, backend)
    return VPolytope(P)
end

"""
    vertices_list(P::HPoly{N};
                  [backend]=default_polyhedra_backend(N),
                  [prunefunc]=removevredundancy!)::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polytope in constraint representation.

### Input

- `P`         -- polytope in constraint representation
- `backend`   -- (optional, default: `default_polyhedra_backend(N)`)
                  the polyhedral computations backend
- `prunefunc` -- (optional, default: `removevredundancy!`) function to post-process
                 the output of `vreps`

### Output

List of vertices.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).

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
function vertices_list(P::HPoly{N};
                       backend=default_polyhedra_backend(N),
                       prunefunc=nothing)::Vector{Vector{N}} where {N<:Real}
    if length(P.constraints) == 0
        return Vector{N}(undef, Vector{N}(undef, 0))
    end
    @assert isdefined(Main, :Polyhedra) "the function `vertices_list` needs " *
                                        "the package 'Polyhedra' to be loaded"
    P = polyhedron(P, backend)
    if prunefunc == nothing
        prunefunc = removevredundancy!
    end
    prunefunc(P)
    return collect(points(P))
end

# ==========================================
# Lower level methods that use Polyhedra.jl
# ==========================================

function load_polyhedra_hpolyhedron() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

function convert(::Type{HPolyhedron{N}}, P::HRep{T, N}) where {T, N}
    constraints = LinearConstraint{N}[]
    for hi in Polyhedra.allhalfspaces(P)
        push!(constraints, HalfSpace(hi.a, hi.β))
    end
    return HPolyhedron(constraints)
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

"""
    polyhedron(P::HPoly{N}, [backend]=default_polyhedra_backend(N)) where {N}

Return an `HRep` polyhedron from `Polyhedra.jl` given a polytope in H-representation.

### Input

- `P`       -- polytope
- `backend` -- (optional, default: call `default_polyhedra_backend(N)`)
                the polyhedral computations backend

### Output

An `HRep` polyhedron.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function polyhedron(P::HPoly{N}, backend=default_polyhedra_backend(N)) where {N}
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
end

end # quote
end # function load_polyhedra_hpolyhedron()
