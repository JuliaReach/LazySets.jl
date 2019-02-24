import Base.rand

export HPolytope,
       vertices_list,
       isbounded

"""
    HPolytope{N<:Real} <: AbstractPolytope{N}

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
struct HPolytope{N<:Real} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N}}

    function HPolytope{N}(constraints::Vector{LinearConstraint{N}};
                          check_boundedness::Bool=false
                         ) where {N<:Real}
        P = new{N}(constraints)
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polytope is not bounded"
        return P
    end
end

# convenience constructor without type parameter
HPolytope(constraints::Vector{LinearConstraint{N}};
          check_boundedness::Bool=false) where {N<:Real} =
    HPolytope{N}(constraints; check_boundedness=check_boundedness)

# constructor with no constraints
HPolytope{N}() where {N<:Real} = HPolytope{N}(Vector{LinearConstraint{N}}())

# constructor with no constraints of type Float64
HPolytope() = HPolytope{Float64}()

# constructor from a simple H-representation
HPolytope(A::AbstractMatrix{N}, b::AbstractVector{N};
          check_boundedness::Bool=false) where {N<:Real} =
    HPolytope(constraints_list(A, b); check_boundedness=check_boundedness)

# constructor from a simple H-representation with type parameter
HPolytope{N}(A::AbstractMatrix{N}, b::AbstractVector{N};
             check_boundedness::Bool=false) where {N<:Real} =
    HPolytope(A, b; check_boundedness=check_boundedness)

promote_rule(::Type{HPolytope{T}},
             ::Type{HPolytope{S}}) where {T<:Real,S<:Real} = HPolytope{promote_type(T,S)}

# --- LazySet interface functions ---


"""
    rand(::Type{HPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::HPolytope{N}

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
              num_vertices::Int=-1
             )::HPolytope{N}
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `rand` needs the " *
                                        "package 'Polyhedra' loaded"
    rng = reseed(rng, seed)
    vpolytope = rand(VPolytope; N=N, dim=dim, rng=rng, seed=seed,
                    num_vertices=num_vertices)
    return convert(HPolytope, vpolytope)
end

"""
    isbounded(P::HPolytope, [use_type_assumption]::Bool=true)::Bool

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
function isbounded(P::HPolytope, use_type_assumption::Bool=true)::Bool
    if use_type_assumption
        return true
    end
    return isbounded(HPolyhedron(P.constraints))
end

function _linear_map_hrep(M::AbstractMatrix{N}, P::HPolytope{N},
                          use_inv::Bool) where {N<:Real}
    constraints = _linear_map_hrep_helper(M, P, use_inv)
    return HPolytope(constraints)
end

# --- functions that use Polyhedra.jl ---


function load_polyhedra_hpolytope() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

@static if VERSION < v"0.7-"

    function convert(::Type{HPolytope{N}}, P::HRep{T, N}) where {T, N}
        constraints = LinearConstraint{N}[]
        for hi in Polyhedra.allhalfspaces(P)
            push!(constraints, HalfSpace(hi.a, hi.β))
        end
        return HPolytope(constraints)
    end

    function convert(::Type{HPolytope}, P::HRep{T, N}) where {T, N}
        return convert(HPolytope{N}, P)
    end

    """
        HPolytope(P::HRep{T, N}) where {T, N}

    Return a polytope in H-representation given a `HRep` polyhedron
    from `Polyhedra.jl`.

    ### Input

    - `P` -- `HRep` polyhedron

    ### Output

    An `HPolytope`.
    """
    function HPolytope(P::HRep{T, N}) where {T, N}
        convert(HPolytope{N}, P)
    end

else

    function convert(::Type{HPolytope{N}}, P::HRep{N}) where {N}
        constraints = LinearConstraint{N}[]
        for hi in Polyhedra.allhalfspaces(P)
            push!(constraints, HalfSpace(hi.a, hi.β))
        end
        return HPolytope(constraints)
    end

    function convert(::Type{HPolytope}, P::HRep{N}) where {N}
        return convert(HPolytope{N}, P)
    end

    function HPolytope(P::HRep{N}) where {N}
        convert(HPolytope{N}, P)
    end

end


end # quote
end # function load_polyhedra_hpolytope()

"""
    vertices_list(P::HPolytope{N};
                  [backend]=default_polyhedra_backend(P, N),
                  [prunefunc]=removevredundancy!)::Vector{Vector{N}} where
                  {N<:Real}

Return the list of vertices of a polytope in constraint representation.

### Input

- `P`         -- polytope in constraint representation
- `backend`   -- (optional, default: `default_polyhedra_backend(P, N)`)
                  the polyhedral computations backend
- `prunefunc` -- (optional, default: `removevredundancy!`) function to
                 post-process the output of `vreps`

### Output

List of vertices.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/).
"""
function vertices_list(P::HPolytope{N};
                       backend=default_polyhedra_backend(P, N),
                       prunefunc=removevredundancy!
                      )::Vector{Vector{N}} where {N<:Real}
    if length(P.constraints) == 0
        return Vector{N}(Vector{N}(undef, 0))
    end
    @assert isdefined(@__MODULE__, :Polyhedra) "the function `vertices_list` needs " *
                                        "the package 'Polyhedra' to be loaded"
    P = polyhedron(P; backend=backend)
    prunefunc(P)
    return collect(points(P))
end

function linear_map(M::Matrix{MN}, P::HPolytope{SN}; kwargs...) where {MN<:Real, SN<:Real}
    T = promote_type(MN, SN)
    MT = convert(Matrix{T}, M)
    PT = convert(HPolytope{T}, P)
    return invoke(linear_map, Tuple{AbstractMatrix{T}, AbstractPolyhedron{T}}, MT, PT; kwargs...)
end
