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
