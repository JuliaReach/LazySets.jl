import Base.rand
import Base.rationalize

export HPolytope,
       vertices_list,
       isbounded

"""
    HPolytope{N, VN<:AbstractVector{N}} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation, that is a finite intersection of half-spaces,

```math
P = \\bigcap_{i = 1}^m H_i,
```
where each ``H_i = \\{x \\in \\mathbb{R}^n : a_i^T x \\leq b_i \\}`` is a half-space,
``a_i \\in \\mathbb{R}^n`` is the normal vector of the ``i``-th half-space and ``b_i`` is the displacement.
It is assumed that ``P`` is bounded (see also [`HPolyhedron`](@ref) which does not make such assumption).

### Fields

- `constraints`       -- vector of linear constraints
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polytope bounded; (boundedness is
                         a running assumption of this type)

### Note

Recall that a polytope is a bounded polyhedron. Boundedness is a running
assumption in this type.
"""
struct HPolytope{N, VN<:AbstractVector{N}} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N, VN}}

    function HPolytope(constraints::Vector{LinearConstraint{N, VN}};
                       check_boundedness::Bool=false) where {N, VN<:AbstractVector{N}}
        P = new{N, VN}(constraints)
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polytope is not bounded"
        return P
    end
end

isoperationtype(::Type{<:HPolytope}) = false
isconvextype(::Type{<:HPolytope}) = true

# constructor for an HPolyhedron with no constraints
function HPolytope{N, VN}() where {N, VN<:AbstractVector{N}}
    HPolytope(Vector{LinearConstraint{N, VN}}())
end

# constructor for an HPolygon with no constraints and given numeric type
function HPolytope{N}() where {N}
    HPolytope(Vector{LinearConstraint{N, Vector{N}}}())
end

# constructor for an HPolytope without explicit numeric type, defaults to Float64
function HPolytope()
    HPolytope{Float64}()
end

# constructor for an HPolytope with constraints of mixed type
function HPolytope(constraints::Vector{<:LinearConstraint})
    HPolytope(_normal_Vector(constraints))
end

# constructor from a simple H-representation
HPolytope(A::AbstractMatrix, b::AbstractVector;
          check_boundedness::Bool=false) =
    HPolytope(constraints_list(A, b); check_boundedness=check_boundedness)


# --- ConvexSet interface functions ---


"""
    rand(::Type{HPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random polytope in constraint representation.

### Input

- `HPolytope`    -- type for dispatch
- `N`            -- (optional, default: `Float64`) numeric type
- `dim`          -- (optional, default: 2) dimension
- `rng`          -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`         -- (optional, default: `nothing`) seed for reseeding
- `num_vertices` -- (optional, default: `-1`) upper bound on the number of
                    vertices of the polytope (see comment below)

### Output

A random polytope in constraint representation.

### Algorithm

We create a random polytope in vertex representation and convert it to
constraint representation (hence the argument `num_vertices`).
See [`rand(::Type{VPolytope})`](@ref).
"""
function rand(::Type{HPolytope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_vertices::Int=-1)
    require(:Polyhedra; fun_name="rand")
    rng = reseed(rng, seed)
    vpolytope = rand(VPolytope; N=N, dim=dim, rng=rng, seed=seed,
                    num_vertices=num_vertices)
    return convert(HPolytope, vpolytope)
end

"""
    isbounded(P::HPolytope, [use_type_assumption]::Bool=true)

Determine whether a polytope in constraint representation is bounded.

### Input

- `P`                   -- polytope in constraint representation
- `use_type_assumption` -- (optional, default: `true`) flag for ignoring the
                           type assumption that polytopes are bounded

### Output

`true` if `use_type_assumption` is activated.
Otherwise, `true` iff `P` is bounded.

### Algorithm

If `!use_type_assumption`, we convert `P` to an `HPolyhedron` `P2` and then use
`isbounded(P2)`.
"""
function isbounded(P::HPolytope, use_type_assumption::Bool=true)
    if use_type_assumption
        return true
    end
    return isbounded(HPolyhedron(P.constraints))
end

function _linear_map_hrep_helper(M::AbstractMatrix, P::HPolytope,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, P, algo)
    return HPolytope(constraints)
end


# --- functions that use Polyhedra.jl ---


function load_polyhedra_hpolytope() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

function convert(::Type{HPolytope}, P::HRep{N}) where {N}
    VT = Polyhedra.hvectortype(P)
    constraints = Vector{LinearConstraint{N, VT}}()
    for hi in Polyhedra.allhalfspaces(P)
        a, b = hi.a, hi.β
        if isapproxzero(norm(a))
            @assert b >= zero(N) "the half-space is inconsistent since it has a " *
                "zero normal direction but the constraint is negative"
            continue
        end
        push!(constraints, HalfSpace(hi.a, hi.β))
    end
    return HPolytope(constraints)
end

"""
    HPolytope(P::HRep)

Return a polytope in H-representation given a `HRep` polyhedron
from `Polyhedra.jl`.

### Input

- `P` -- `HRep` polyhedron

### Output

An `HPolytope`.
"""
function HPolytope(P::HRep)
    convert(HPolytope, P)
end

end # quote
end # function load_polyhedra_hpolytope()

"""
    vertices_list(P::HPolytope{N};
                  [backend]=nothing, [prune]::Bool=true) where {N}

Return the list of vertices of a polytope in constraint representation.

### Input

- `P`       -- polytope in constraint representation
- `backend` -- (optional, default: `nothing`) the polyhedral computations backend
- `prune`   -- (optional, default: `true`) flag to remove redundant vertices

### Output

List of vertices.

### Algorithm

If the polytope is two-dimensional, the polytope is converted to a polygon in
H-representation and then its `vertices_list` function is used. This ensures
that, by default, the optimized two-dimensional methods are used.

It is possible to use the `Polyhedra` backend in two-dimensions as well
by passing, e.g. `backend=CDDLib.Library()`.

If the polytope is not two-dimensional, the concrete polyhedra manipulation
library `Polyhedra` is used. The actual computation is performed by a given
backend; for the default backend used in `LazySets` see `default_polyhedra_backend(P)`.
For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function vertices_list(P::HPolytope{N};
                       backend=nothing, prune::Bool=true) where {N}
    if length(P.constraints) == 0
        return Vector{N}(Vector{N}(undef, 0))
    end

    if dim(P) == 2 && backend == nothing
        return vertices_list(convert(HPolygon, P, prune=prune))
    else
        require(:Polyhedra; fun_name="vertices_list")
        if backend == nothing
            backend = default_polyhedra_backend(P)
        end
        Q = polyhedron(P; backend=backend)
        if prune
            removevredundancy!(Q; ztol=_ztol(N))
        end
        return collect(Polyhedra.points(Q))
    end
end

# used for dispatch, see minkowski_sum(::AbstractPolytope{N}, ::AbstractPolytope{N}; ...)
function _vertices_list(P::HPolytope, backend)
    return vertices_list(P, backend=backend)
end

# ============================================
# Functionality that requires Symbolics
# ============================================
function load_symbolics_hpolytope()

return quote

"""
    HPolytope(expr::Vector{<:Num}, vars=_get_variables(expr); [N]::Type{<:Real}=Float64, [check_boundedness]::Bool=false)

Return the polytope in half-space representation given by a list of symbolic expressions.

### Input

- `expr` -- vector of symbolic expressions that describes each half-space
- `vars` -- (optional, default: `_get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may be incorrect; it is advised
            to always pass `vars` explicitly)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned half-space
- `check_boundedness` -- (optional, default: `false`) flag to check boundedness

### Output

An `HPolytope`.

### Examples

```jldoctest
julia> using Symbolics

julia> vars = @variables x y
2-element Vector{Num}:
 x
 y

julia> HPolytope([x <= 1, x >= 0, y <= 1, y >= 0], vars)
HPolytope{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, 0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([0.0, -1.0], 0.0)])
```
"""
function HPolytope(expr::Vector{<:Num}, vars::AbstractVector{Num};
                   N::Type{<:Real}=Float64, check_boundedness::Bool=false)
    return HPolytope([HalfSpace(ex, vars; N=N) for ex in expr],
                     check_boundedness=check_boundedness)
end

HPolytope(expr::Vector{<:Num}; N::Type{<:Real}=Float64, check_boundedness::Bool=false) = HPolytope(expr, _get_variables(expr); N=N, check_boundedness=check_boundedness)
HPolytope(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64, check_boundedness::Bool=false) = HPolytope(expr, _vec(vars); N=N, check_boundedness=check_boundedness)

end end  # quote / load_modeling_toolkit_hpolytope()
