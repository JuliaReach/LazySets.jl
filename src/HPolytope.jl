import Base.rand

export HPolytope,
       vertices_list

"""
    HPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation.

### Fields

- `constraints` -- vector of linear constraints

### Note

Recall that a polytope is a bounded polyhedron. Boundedness is a running
assumption in this type.
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
function HPolytope(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[]
    @inbounds for i in 1:m
        push!(constraints, LinearConstraint(A[i, :], b[i]))
    end
    return HPolytope(constraints)
end

HPolytope{N}(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real} = HPolytope(A, b)

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
