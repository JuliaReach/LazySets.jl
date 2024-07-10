export Universe

"""
    Universe{N} <: AbstractPolyhedron{N}

Type that represents the universal set, i.e., the set of all elements.

### Fields

- `dim` -- the ambient dimension of the set
"""
struct Universe{N} <: AbstractPolyhedron{N}
    dim::Int
end

isoperationtype(::Type{<:Universe}) = false

# default constructor of type Float64
Universe(dim::Int) = Universe{Float64}(dim)

"""
    constraints(U::Universe{N}) where {N}

Construct an iterator over the constraints of a universe.

### Input

- `U` -- universe

### Output

The empty iterator, as the universe is unconstrained.
"""
function constraints(U::Universe{N}) where {N}
    return EmptyIterator{Vector{N}}()
end

"""
    constraints_list(U::Universe{N}) where {N}

Return the list of constraints defining a universe.

### Input

- `U` -- universe

### Output

The empty list of constraints, as the universe is unconstrained.
"""
function constraints_list(U::Universe{N}) where {N}
    return HalfSpace{N,Vector{N}}[]
end

"""
    constrained_dimensions(U::Universe)

Return the indices in which a universe is constrained.

### Input

- `U` -- universe

### Output

The empty vector, as the universe is unconstrained in every dimension.
"""
function constrained_dimensions(U::Universe)
    return Int[]
end

"""
    dim(U::Universe)

Return the dimension of a universe.

### Input

- `U` -- universe

### Output

The ambient dimension of a universe.
"""
function dim(U::Universe)
    return U.dim
end

"""
    ρ(d::AbstractVector, U::Universe)

Return the support function of a universe.

### Input

- `d` -- direction
- `U` -- universe

### Output

The support function in the given direction.

### Algorithm

If the direction is all zero, the result is zero.
Otherwise, the result is `Inf`.
"""
function ρ(d::AbstractVector, U::Universe)
    N = promote_type(eltype(d), eltype(U))
    return iszero(d) ? zero(N) : N(Inf)
end

"""
    σ(d::AbstractVector, U::Universe)

Return the support vector of a universe.

### Input

- `d` -- direction
- `U` -- universe

### Output

A vector with infinity values, except in dimensions where the direction is zero.
"""
function σ(d::AbstractVector, U::Universe)
    N = promote_type(eltype(d), eltype(U))
    return [iszero(v) ? v : v > zero(N) ? N(Inf) : N(-Inf) for v in d]
end

"""
    ∈(x::AbstractVector, U::Universe)

Check whether a given point is contained in a universe.

### Input

- `x` -- point/vector
- `U` -- universe

### Output

`true`.

### Examples

```jldoctest
julia> [1.0, 0.0] ∈ Universe(2)
true
```
"""
function ∈(x::AbstractVector, U::Universe)
    @assert length(x) == dim(U) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(U))-dimensional set"
    return true
end

"""
    an_element(U::Universe{N}) where {N}

Return some element of a universe.

### Input

- `U` -- universe

### Output

The origin.
"""
function an_element(U::Universe{N}) where {N}
    return zeros(N, dim(U))
end

"""
    rand(::Type{Universe}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a universe (note that there is nothing to randomize).

### Input

- `Universe` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 2) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

The (only) universe of the given numeric type and dimension.
"""
function rand(::Type{Universe};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return Universe{N}(dim)
end

"""
    isempty(U::Universe)

Check whether a universe is empty.

### Input

- `U` -- universe

### Output

`false`.
"""
function isempty(U::Universe)
    return false
end

"""
    isbounded(U::Universe)

Check whether a universe is bounded.

### Input

- `U` -- universe

### Output

`false`.
"""
function isbounded(U::Universe)
    return false
end

function isboundedtype(::Type{<:Universe})
    return false
end

"""
    isuniversal(U::Universe{N}, [witness]::Bool=false) where {N}

Check whether a universe is universal.

### Input

- `U`       -- universe
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true`
* If `witness` option is activated: `(true, [])`
"""
function isuniversal(U::Universe{N}, witness::Bool=false) where {N}
    return witness ? (true, N[]) : true
end

"""
    norm(U::Universe, [p]::Real=Inf)

Return the norm of a universe.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(U::Universe, p::Real=Inf)
    return error("a universe does not have a norm")
end

"""
    radius(U::Universe, [p]::Real=Inf)

Return the radius of a universe.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(U::Universe, p::Real=Inf)
    return error("a universe does not have a radius")
end

"""
    diameter(U::Universe, [p]::Real=Inf)

Return the diameter of a universe.
It is the diameter of the enclosing ball (of the given ``p``-norm) of minimal
volume .

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(U::Universe, p::Real=Inf)
    return error("a universe does not have a diameter")
end

"""
    translate(U::Universe, v::AbstractVector)

Translate (i.e., shift) a universe by a given vector.

### Input

- `U` -- universe
- `v` -- translation vector

### Output

The universe.
"""
function translate(U::Universe, v::AbstractVector)
    return translate!(U, v)  # no need to copy
end

"""
    translate!(U::Universe, v::AbstractVector)

Translate (i.e., shift) a universe by a given vector, in-place.

### Input

- `U` -- universe
- `v` -- translation vector

### Output

The universe.
"""
function translate!(U::Universe, v::AbstractVector)
    @assert length(v) == dim(U) "cannot translate a $(dim(U))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return U
end

function linear_map_inverse(Minv::AbstractMatrix{N}, U::Universe{N}) where {N}
    @assert size(Minv, 1) == dim(U) "a linear map of size $(size(Minv)) " *
                                    "cannot be applied to a universe of dimension $(dim(U))"
    n = size(Minv, 2)
    return Universe{N}(n)
end

function project(U::Universe{N}, block::AbstractVector{Int}; kwargs...) where {N}
    return Universe{N}(length(block))
end

"""
    permute(U::Universe, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `U` -- universe
- `p` -- permutation vector

### Output

The same universe.
"""
function permute(U::Universe, p::AbstractVector{Int})
    return U
end

function tosimplehrep(U::Universe)
    return tosimplehrep(constraints_list(U); n=dim(U))
end

"""
    complement(U::Universe{N}) where {N}

Return the complement of an universe.

### Input

- `∅` -- universe

### Output

The empty set of the same dimension.
"""
function complement(U::Universe{N}) where {N}
    return EmptySet{N}(dim(U))
end

function load_polyhedra_universe() # function to be loaded by Requires
    return quote
        # see the interface file init_Polyhedra.jl for the imports

        """
            polyhedron(U::Universe; [backend]=default_polyhedra_backend(P))

        Return an `HRep` polyhedron from `Polyhedra.jl` given a universe.

        ### Input

        - `U`       -- universe
        - `backend` -- (optional, default: call `default_polyhedra_backend(P)`)
                       the backend for polyhedral computations

        ### Output

        An `HRep` polyhedron.

        ### Notes

        For further information on the supported backends see
        [Polyhedra's documentation](https://juliapolyhedra.github.io/).
        """
        function polyhedron(U::Universe; backend=default_polyhedra_backend(U))
            A, b = tosimplehrep(U)
            return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
        end
    end
end  # quote / load_polyhedra_universe()

"""
    reflect(U::Universe)

Concrete reflection of a universe `U`, resulting in the reflected set `-U`.

### Input

- `U` -- universe

### Output

The same universe.
"""
function reflect(U::Universe)
    return U
end

function scale(α::Real, U::Universe{N}) where {N}
    if iszero(α)
        return ZeroSet{N}(dim(U))
    end
    return U
end

function scale!(α::Real, U::Universe)
    if iszero(α)
        throw(ArgumentError("cannot 0-scale a universe in-place"))
    end
    return U
end
