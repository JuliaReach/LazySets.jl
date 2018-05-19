import LazySets.dim

"""
    AbstractDirections{N}

Abstract type for template direction representations.

### Notes

All subtypes should implement the standard iterator methods from `Base` and the
function `dim(d<:AbstractDirections)::Int`.
"""
abstract type AbstractDirections{N} end

"""
    dim(ad::AbstractDirections)::Int

Returns the dimension of the generated directions.

### Input

- `ad` -- box direction representation

### Output

The dimension of the generated directions.
"""
function dim(ad::AbstractDirections)::Int
    return ad.n
end

# box directions

"""
    UnitVector{T} <: AbstractVector{T}

A lazy unit vector with arbitrary one-element.

### Fields

- `i` -- index of non-zero entry
- `n` -- vector length
- `v` -- non-zero entry
"""
struct UnitVector{T} <: AbstractVector{T}
    i::Int
    n::Int
    v::T
end

function Base.getindex(e::UnitVector{T}, i::Int) where T
    @boundscheck @assert 1 <= i <= e.n
    return i == e.i ? e.v : zero(T)
end

Base.size(e::UnitVector) = (e.n,)

"""
    BoxDirections{N} <: AbstractDirections{N}

Box direction representation.

### Fields

- `n` -- dimension
"""
struct BoxDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
BoxDirections(n::Int) = BoxDirections{Float64}(n)

Base.eltype(::Type{BoxDirections{N}}) where N = AbstractVector{N}
Base.start(bd::BoxDirections) = 1
Base.next(bd::BoxDirections{N}, state) where N = (
    UnitVector{N}(abs(state), bd.n, convert(N, sign(state))), # value
    state == bd.n ? -bd.n : state + 1) # next state
Base.done(bd::BoxDirections, state) = state == 0
Base.length(bd::BoxDirections) = 2*bd.n

# box-diagonal directions

"""
    BoxDiagDirections{N} <: AbstractDirections{N}

Box-diagonal direction representation.

### Fields

- `n` -- dimension

### Notes

Box-diagonal directions can be seen as the union of diagonal directions (all
entries are ±1) and box directions (one entry is ±1, all other entries are 0).
The iterator first enumerates all diagonal directions, and then all box
directions.
"""
struct BoxDiagDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
BoxDiagDirections(n::Int) = BoxDiagDirections{Float64}(n)

Base.eltype(::Type{BoxDiagDirections{N}}) where N = AbstractVector{N}
Base.start(bdd::BoxDiagDirections{N}) where N = ones(N, bdd.n)
function Base.next(bdd::BoxDiagDirections{N}, state::AbstractVector) where N
    i = 1
    while i <= bdd.n && state[i] < 0
        state[i] = -state[i]
        i = i+1
    end
    if i > bdd.n
        if bdd.n == 1
            return (copy(state), 0) # finish here to avoid duplicates
        else
            return (copy(state), 1) # continue with box directions
        end
    else
        state[i] = -state[i]
        return (copy(state), state)
    end
end
function Base.next(bdd::BoxDiagDirections{N}, state::Int) where N
    return next(BoxDirections{N}(bdd.n), state)
end
Base.done(bdd::BoxDiagDirections, state) = state == 0
Base.length(bdd::BoxDiagDirections) = bdd.n == 1 ? 2 : 2^bdd.n + 2*bdd.n
