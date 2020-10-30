import Base.rand

export HPolytope,
       vertices_list,
       isbounded

"""
    HPolytope{N<:Real, VN<:AbstractVector{N}} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation.

### Fields

- `constraints`       -- vector of linear constraints
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polytope bounded; (boundedness is
                         a running assumption of this type)

### Note

Recall that a polytope is a bounded polyhedron. Boundedness is a running
assumption in this type.
"""
struct HPolytope{N<:Real, VN<:AbstractVector{N}} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N, VN}}

    function HPolytope(constraints::Vector{LinearConstraint{N, VN}};
                       check_boundedness::Bool=false) where {N<:Real, VN<:AbstractVector{N}}
        P = new{N, VN}(constraints)
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polytope is not bounded"
        return P
    end
end

isoperationtype(::Type{<:HPolytope}) = false
isconvextype(::Type{<:HPolytope}) = true

# constructor for an HPolyhedron with no constraints
function HPolytope{N, VN}() where {N<:Real, VN<:AbstractVector{N}}
    HPolytope(Vector{LinearConstraint{N, VN}}())
end

# constructor for an HPolygon with no constraints and given numeric type
function HPolytope{N}() where {N<:Real}
    HPolytope(Vector{LinearConstraint{N, Vector{N}}}())
end

# constructor for an HPolytope without explicit numeric type, defaults to Float64
function HPolytope()
    HPolytope{Float64}()
end

# constructor from a simple H-representation
HPolytope(A::AbstractMatrix{N}, b::AbstractVector{N};
          check_boundedness::Bool=false) where {N<:Real} =
    HPolytope(constraints_list(A, b); check_boundedness=check_boundedness)


# --- LazySet interface functions ---


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

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::HPolytope{N},
                                 algo::AbstractLinearMapAlgorithm) where {N<:Real}
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
    HPolytope(P::HRep{N}) where {N}

Return a polytope in H-representation given a `HRep` polyhedron
from `Polyhedra.jl`.

### Input

- `P` -- `HRep` polyhedron

### Output

An `HPolytope`.
"""
function HPolytope(P::HRep{N}) where {N}
    convert(HPolytope, P)
end

end # quote
end # function load_polyhedra_hpolytope()

"""
    vertices_list(P::HPolytope{N};
                  [backend]=nothing, [prune]::Bool=true) where {N<:Real}

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
backend; for the default backend used in `LazySets` see `default_polyhedra_backend(N)`.
For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function vertices_list(P::HPolytope{N};
                       backend=nothing, prune::Bool=true) where {N<:Real}
    if length(P.constraints) == 0
        return Vector{N}(Vector{N}(undef, 0))
    end

    if dim(P) == 2 && backend == nothing
        return vertices_list(convert(HPolygon, P, prune=prune))
    else
        require(:Polyhedra; fun_name="vertices_list")
        if backend == nothing
            backend = default_polyhedra_backend(P, N)
        end
        Q = polyhedron(P; backend=backend)
        if prune
            removevredundancy!(Q)
        end
        return collect(Polyhedra.points(Q))
    end
end

# used for dispatch, see minkowski_sum(::AbstractPolytope{N}, ::AbstractPolytope{N}; ...)
function _vertices_list(P::HPolytope, backend)
    return vertices_list(P, backend=backend)
end

# ============================================
# Functionality that requires ModelingToolkit
# ============================================
function load_modeling_toolkit_hpolytope()

return quote

"""
    HPolytope(expr::Vector{<:Operation}, vars=get_variables(first(expr)); N::Type{<:Real}=Float64)

Return the polytope in half-space representation given by a list of symbolic expressions.

### Input

- `expr` -- vector of symbolic expressions that describes each half-space
- `vars` -- (optional, default: `get_variables(expr)`), if an array of variables is given,
            use those as the ambient variables in the set with respect to which derivations
            take place; otherwise, use only the variables which appear in the given
            expression (but be careful because the order may change)
- `N`    -- (optional, default: `Float64`) the numeric type of the returned half-space

### Output

An `HPolytope`.

### Examples

```julia
julia> using ModelingToolkit

julia> vars = @variables x y
(x, y)

julia> HPolytope([x <= 1, x >= 0, y <= 1, y >= 0], vars)
HPolytope{Float64,Array{Float64,1}}(HalfSpace{Float64,Array{Float64,1}}[HalfSpace{Float64,Array{Float64,1}}([1.0, 0.0], 1.0),
HalfSpace{Float64,Array{Float64,1}}([-1.0, 0.0], 0.0), HalfSpace{Float64,Array{Float64,1}}([0.0, 1.0], 1.0),
HalfSpace{Float64,Array{Float64,1}}([0.0, -1.0], 0.0)])
```
"""
function HPolytope(expr::Vector{<:Operation}, vars=get_variables(first(expr));
                   N::Type{<:Real}=Float64, check_boundedness::Bool=false)
    return HPolytope([HalfSpace(ex, vars; N=N) for ex in expr], check_boundedness=check_boundedness)
end

end end  # quote / load_modeling_toolkit_hpolytope()
