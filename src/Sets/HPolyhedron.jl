import Base: isempty,
             rand,
             convert

export HPolyhedron,
       dim, σ, ∈,
       addconstraint!,
       constraints_list,
       tohrep, tovrep,
       convex_hull,
       cartesian_product,
       isempty,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       constrained_dimensions

"""
    HPolyhedron{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a convex polyhedron in H-representation.

### Fields

- `constraints` -- vector of linear constraints
"""
struct HPolyhedron{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    constraints::Vector{LinearConstraint{N, VN}}

    function HPolyhedron(constraints::Vector{LinearConstraint{N, VN}}) where {N<:Real,
                                                                              VN<:AbstractVector{N}}
        return new{N, VN}(constraints)
    end
end

isoperationtype(::Type{<:HPolyhedron}) = false
isconvextype(::Type{<:HPolyhedron}) = true

# constructor for an HPolyhedron with no constraints
function HPolyhedron{N, VN}() where {N<:Real, VN<:AbstractVector{N}}
    HPolyhedron(Vector{LinearConstraint{N, VN}}())
end

# constructor for an HPolyhedron with no constraints and given numeric type
function HPolyhedron{N}() where {N<:Real}
    HPolyhedron(Vector{LinearConstraint{N, Vector{N}}}())
end

# constructor for an HPolyhedron without explicit numeric type, defaults to Float64
function HPolyhedron()
    HPolyhedron{Float64}()
end

# constructor from a simple H-representation
HPolyhedron(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real} =
    HPolyhedron(constraints_list(A, b))

# convenience union type
const HPoly{N} = Union{HPolytope{N}, HPolyhedron{N}}


# --- LazySet interface functions ---


"""
    dim(P::HPoly{N}) where {N<:Real}

Return the dimension of a polyhedron in H-representation.

### Input

- `P`  -- polyhedron in H-representation

### Output

The ambient dimension of the polyhedron in H-representation.
If it has no constraints, the result is ``-1``.
"""
function dim(P::HPoly{N}) where {N<:Real}
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    ρ(d::AbstractVector{N}, P::HPoly{N};
      solver=default_lp_solver(N)) where {N<:Real}

Evaluate the support function of a polyhedron (in H-representation) in a given
direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in H-representation
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Output

The support function of the polyhedron.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, P::HPoly{N};
           solver=default_lp_solver(N)) where {N<:Real}
    lp, unbounded = σ_helper(d, P, solver)
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
    σ(d::AbstractVector{N}, P::HPoly{N}; solver=default_lp_solver(N)
     ) where {N<:Real}

Return the support vector of a polyhedron (in H-representation) in a given
direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in H-representation
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, P::HPoly{N}; solver=default_lp_solver(N)
          ) where {N<:Real}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            error("the support vector in direction $(d) is undefined because " *
                  "the polytope is unbounded")
        end
        # construct the solution from the solver's ray result
        if lp == nothing
            ray = d
        elseif haskey(lp.attrs, :unboundedray)
            ray = lp.attrs[:unboundedray]
        else
            error("LP solver did not return an infeasibility ray")
        end
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

function σ_helper(d::AbstractVector{N}, P::HPoly{N}, solver) where {N<:Real}
    # let c = -d as a Vector since GLPK does not accept sparse vectors
    # (see #1011)
    c = to_negative_vector(d)

    (A, b) = tosimplehrep(P)
    if length(b) == 0
        unbounded = true
        lp = nothing
    else
        sense = '<'
        l = -Inf
        u = Inf
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
    rand(::Type{HPolyhedron}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    P = rand(HPolytope; N=N, dim=dim, rng=rng)
    constraints_P = constraints_list(P)
    constraints_Q = Vector{eltype(constraints_P)}()
    for i in 1:length(constraints_P)
        if rand(Bool)
            push!(constraints_Q, constraints_P[i])
        end
    end
    return HPolyhedron(constraints_Q)
end


# ===========================================
# HPolyhedron and HPolytope's shared methods
# ===========================================


"""
    addconstraint!(P::HPoly{N}, constraint::LinearConstraint{N}) where {N<:Real}

Add a linear constraint to a polyhedron in H-representation.

### Input

- `P`          -- polyhedron in H-representation
- `constraint` -- linear constraint to add

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(P::HPoly{N},
                        constraint::LinearConstraint{N}) where {N<:Real}
    push!(P.constraints, constraint)
    return nothing
end

"""
    constraints_list(P::HPoly{N}) where {N<:Real}

Return the list of constraints defining a polyhedron in H-representation.

### Input

- `P` -- polyhedron in H-representation

### Output

The list of constraints of the polyhedron.
"""
function constraints_list(P::HPoly{N}) where {N<:Real}
    return P.constraints
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
    normalize(P::HPoly{N}, p=N(2)) where {N<:Real}

Normalize a polyhedron in constraint representation.

### Input

- `P` -- polyhedron in constraint representation
- `p` -- (optional, default: `2`) norm

### Output

A new polyhedron in constraint representation whose normal directions ``a_i``
are normalized, i.e., such that ``‖a_i‖_p = 1`` holds.
"""
function normalize(P::HPoly{N}, p=N(2)) where {N<:Real}
    constraints = [normalize(hs, p) for hs in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

"""
    remove_redundant_constraints(P::HPoly{N};
                                 backend=default_lp_solver(N)
                                ) where {N<:Real}

Remove the redundant constraints in a polyhedron in H-representation.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `default_lp_solver(N)`) the backend used to
               solve the linear program

### Output

A polyhedron equivalent to `P` but with no redundant constraints, or an empty
set if `P` is detected to be empty, which may happen if the constraints are
infeasible.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:LinearConstraint})` for
details.
"""
function remove_redundant_constraints(P::HPoly{N};
                                      backend=default_lp_solver(N)
                                     ) where {N<:Real}
    Pred = copy(P)
    if remove_redundant_constraints!(Pred, backend=backend)
        return Pred
    else # the polyhedron P is empty
        return EmptySet{N}(dim(P))
    end
end

"""
    remove_redundant_constraints!(P::HPoly{N};
                                  backend=default_lp_solver(N)) where {N<:Real}

Remove the redundant constraints in a polyhedron in H-representation; the
polyhedron is updated in-place.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `default_lp_solver(N)`) the backend used to
               solve the linear program

### Output

`true` if the method was successful and the polyhedron `P` is modified by
removing its redundant constraints, and `false` if `P` is detected to be empty,
which may happen if the constraints are infeasible.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:LinearConstraint})` for
details.
"""
function remove_redundant_constraints!(P::HPoly{N};
                                       backend=default_lp_solver(N)
                                      ) where {N<:Real}
    remove_redundant_constraints!(P.constraints, backend=backend)
end

"""
    translate(P::HPoly{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a polyhedron in constraint representation by a given
vector.

### Input

- `P`     -- polyhedron in constraint representation
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated polyhedron in constraint representation.

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.

### Algorithm

We translate every constraint.
"""
function translate(P::HPoly{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

# ========================================================
# External methods that require Polyhedra.jl to be loaded
# ========================================================

"""
    convex_hull(P1::HPoly{N}, P2::HPoly{N};
               [backend]=default_polyhedra_backend(P1, N)) where {N<:Real}

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
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function convex_hull(P1::HPoly{N},
                     P2::HPoly{N};
                     backend=default_polyhedra_backend(P1, N)) where {N<:Real}
    require(:Polyhedra; fun_name="convex_hull")
    Pch = Polyhedra.convexhull(polyhedron(P1; backend=backend),
                               polyhedron(P2; backend=backend))
    removehredundancy!(Pch)
    return convert(basetype(P1), Pch)
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
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function cartesian_product(P1::HPoly{N},
                           P2::HPoly{N};
                           backend=default_polyhedra_backend(P1, N)
                          ) where {N<:Real}
    require(:Polyhedra; fun_name="`cartesian_product")
    Pcp = Polyhedra.hcartesianproduct(polyhedron(P1; backend=backend),
                                      polyhedron(P2; backend=backend))
    return convert(basetype(P1), Pcp)
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
For further information on the supported backends see [Polyhedra's
documentation](https://juliapolyhedra.github.io/).
"""
function tovrep(P::HPoly{N};
                backend=default_polyhedra_backend(P, N)) where {N<:Real}
    require(:Polyhedra; fun_name="tovrep")
    P = polyhedron(P; backend=backend)
    return VPolytope(P)
end

"""
   isempty(P::HPoly{N}, witness::Bool=false;
           [use_polyhedra_interface]::Bool=false, [solver]=default_lp_solver(N),
           [backend]=nothing) where {N<:Real}

Determine whether a polyhedron is empty.

### Input

- `P`       -- polyhedron
- `witness` -- (optional, default: `false`) compute a witness if activated
- `use_polyhedra_interface` -- (optional, default: `false`) if `true`, we use
               the `Polyhedra` interface for the emptiness test
- `solver`  -- (optional, default: `default_lp_solver(N)`) LP-solver backend
- `backend` -- (optional, default: `nothing`) backend for polyhedral
               computations in `Polyhedra`; its value is set internally (see the
               Notes below for details)

### Output

* If `witness` option is deactivated: `true` iff ``P = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``P = ∅``
  * `(false, v)` iff ``P ≠ ∅`` and ``v ∈ P``

### Notes

The default value of the `backend` is set internally and depends on whether the
`use_polyhedra_interface` option is set or not.
If the option is set, we use `default_polyhedra_backend(P, N)`.

Witness production is not supported if `use_polyhedra_interface` is `true`.

### Algorithm

The algorithm sets up a feasibility LP for the constraints of `P`.
If `use_polyhedra_interface` is `true`, we call `Polyhedra.isempty`.
Otherwise, we set up the LP internally.
"""
function isempty(P::HPoly{N},
                 witness::Bool=false;
                 use_polyhedra_interface::Bool=false,
                 solver=default_lp_solver(N),
                 backend=nothing) where {N<:Real}
    if length(constraints_list(P)) < 2
        # catch corner case because of problems in LP solver for Rationals
        return witness ? (false, an_element(P)) : false
    end
    if use_polyhedra_interface
        require(:Polyhedra; fun_name="isempty", explanation="with the active " *
            "option `use_polyhedra_interface`")
        if backend == nothing
            backend = default_polyhedra_backend(P, N)
        end
        result = Polyhedra.isempty(polyhedron(P; backend=backend), solver)
        if result
            return witness ? (true, N[]) : true
        elseif witness
            error("witness production is not supported yet")
        else
            return false
        end
    else
        A, b = tosimplehrep(P)
        lbounds, ubounds = -Inf, Inf
        sense = '<'
        obj = zeros(N, size(A, 2))
        lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
        if lp.status == :Optimal
            return witness ? (false, lp.sol) : false
        elseif lp.status == :Infeasible
            return witness ? (true, N[]) : true
        end
        error("LP returned status $(lp.status) unexpectedly")
    end
end

convert(::Type{HPolytope}, P::HPolyhedron{N}) where {N<:Real} =
    HPolytope(copy(constraints_list(P)))
convert(::Type{HPolyhedron}, P::HPolytope{N}) where {N<:Real} =
    HPolyhedron(copy(constraints_list(P)))

# ==========================================
# Lower level methods that use Polyhedra.jl
# ==========================================

function load_polyhedra_hpolyhedron() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

function convert(::Type{HPolyhedron}, P::HRep{N}) where {N}
    VN = Polyhedra.hvectortype(P)
    constraints = Vector{LinearConstraint{N, VN}}()
    for hi in Polyhedra.allhalfspaces(P)
        a, b = hi.a, hi.β
        if isapproxzero(norm(a))
            continue
        end
        push!(constraints, HalfSpace(a, b))
    end
    return HPolyhedron(constraints)
end

"""
     HPolyhedron(P::HRep{N}) where {N}

Return a polyhedron in H-representation given a `HRep` polyhedron
from `Polyhedra.jl`.

### Input

- `P` -- `HRep` polyhedron

### Output

An `HPolyhedron`.
"""
function HPolyhedron(P::HRep{N}) where {N}
    convert(HPolyhedron, P)
end

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
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function polyhedron(P::HPoly{N};
                    backend=default_polyhedra_backend(P, N)) where {N<:Real}
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
end

end # quote
end # function load_polyhedra_hpolyhedron()

"""
    _isbounded_stiemke(P::HPolyhedron{N}; solver=LazySets.default_lp_solver(N),
                       check_nonempty::Bool=true) where {N<:Real}

Determine whether a polyhedron is bounded using Stiemke's theorem of alternatives.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `default_lp_solver(N)`) the backend used
               to solve the linear program
- `check_nonempty` -- (optional, default: `true`) if `true`, check the
                      precondition to this algorithm that `P` is non-empty

### Output

`true` iff the polyhedron is bounded

### Notes

The algorithm internally calls `isempty` to check whether the polyhedron is empty.
This computation can be avoided using the `check_nonempty` flag.

### Algorithm

The algorithm is based on Stiemke's theorem of alternatives, see e.g. [1].

Let the polyhedron ``P`` be given in constraint form ``Ax ≤ b``. We assume that
the polyhedron is not empty.

Proposition 1. If ``\\ker(A)≠\\{0\\}``, then ``P`` is unbounded.

Proposition 2. Assume that ``ker(A)={0}`` and ``P`` is non-empty.
Then ``P`` is bounded if and only if the following linear
program admits a feasible solution: ``\\min∥y∥_1`` subject to ``A^Ty=0`` and ``y≥1``.

[1] Mangasarian, Olvi L. *Nonlinear programming.*
    Society for Industrial and Applied Mathematics, 1994.
"""
function _isbounded_stiemke(P::HPolyhedron{N}; solver=LazySets.default_lp_solver(N),
                            check_nonempty::Bool=true) where {N<:Real}
    if check_nonempty && isempty(P)
        return true
    end

    A, b = tosimplehrep(P)
    m, n = size(A)

    if !isempty(nullspace(A))
        return false
    end

    At = copy(transpose(A))
    c = ones(N, m)
    lp = linprog(c, At, '=', zeros(n), one(N), Inf, solver)
    return (lp.status == :Optimal)
end

# ============================================
# Functionality that requires ModelingToolkit
# ============================================
function load_modeling_toolkit_hpolyhedron()

return quote

"""
    HPolyhedron(expr::Vector{<:Operation}, vars=get_variables(first(expr)); N::Type{<:Real}=Float64)

Return the polyhedron in half-space representation given by a list of symbolic expressions.

### Input

- `expr` -- vector of symbolic expressions that describes each half-space
- `vars` -- (optional, default: `get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may change)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned half-space

### Output

An `HPolyhedron`.

### Examples

```julia
julia> using ModelingToolkit

julia> vars = @variables x y
(x, y)

julia> HPolyhedron([x + y <= 1, x + y >= -1], vars)
HPolyhedron{Float64,Array{Float64,1}}(HalfSpace{Float64,Array{Float64,1}}[HalfSpace{Float64,Array{Float64,1}}([1.0, 1.0], 1.
0), HalfSpace{Float64,Array{Float64,1}}([-1.0, -1.0], 1.0)])
```
"""
function HPolyhedron(expr::Vector{<:Operation}, vars=get_variables(first(expr)); N::Type{<:Real}=Float64)
    return HPolyhedron([HalfSpace(ex, vars; N=N) for ex in expr])
end

end end  # quote / load_modeling_toolkit_hpolyhedron()
