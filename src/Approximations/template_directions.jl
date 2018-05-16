import LazySets.dim

"""
    AbstractDirections{N}

Abstract type for template direction representations.

### Notes

All subtypes should implement the standard iterator methods from `Base` and the
function `dim(d<:AbstractDirections)::Int`.
"""
abstract type AbstractDirections{N} end

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

"""
    dim(bd::BoxDirections)::Int

Returns the dimension of the generated directions.

### Input

- `bd` -- box direction representation

### Output

The dimension of the generated directions.
"""
function dim(bd::BoxDirections)::Int
    return bd.n
end

# octagon directions

"""
    OctDirections{N} <: AbstractDirections{N}

Octagon direction representation.

### Fields

- `n` -- dimension

### Notes

Octagon directions can be seen as the union of diagonal directions (all entries
are ±1) and box directions (one entry is ±1, all other entries are 0).
The iterator first enumerates all diagonal directions; then it enumerates all
box directions.
"""
struct OctDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
OctDirections(n::Int) = OctDirections{Float64}(n)

Base.eltype(::Type{OctDirections{N}}) where N = AbstractVector{N}
Base.start(od::OctDirections{N}) where N = ones(N, od.n)
function Base.next(od::OctDirections{N}, state::AbstractVector) where N
    i = 1
    while i <= od.n && state[i] < 0
        state[i] = -state[i]
        i = i+1
    end
    if i > od.n
        if od.n == 1
            return (copy(state), 0) # finish here to avoid duplicates
        else
            return (copy(state), 1) # continue with box directions
        end
    else
        state[i] = -state[i]
        return (copy(state), state)
    end
end
function Base.next(od::OctDirections{N}, state::Int) where N
    return next(BoxDirections{N}(od.n), state)
end
Base.done(od::OctDirections, state) = state == 0
Base.length(od::OctDirections) = od.n == 1 ? 2 : 2^od.n + 2*od.n

"""
    dim(od::OctDirections)::Int

Returns the dimension of the generated directions.

### Input

- `od` -- octagon direction representation

### Output

The dimension of the generated directions.
"""
function dim(od::OctDirections)::Int
    return od.n
end
