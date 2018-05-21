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
Base.length(bd::BoxDirections) = 2 * bd.n

# octagon directions

"""
    OctDirections{N} <: AbstractDirections{N}

Octagon direction representation.

### Fields

- `n` -- dimension

### Notes

Octagon directions consist of all vectors that are zero almost everywhere except
in two dimensions ``i``, ``j`` (possibly ``i = j``) where it is ``±1``.
"""
struct OctDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
OctDirections(n::Int) = OctDirections{Float64}(n)

Base.eltype(::Type{OctDirections{N}}) where N = AbstractVector{N}
function Base.start(od::OctDirections{N}) where N
    if od.n == 1
        return (1)
    end
    vec = zeros(N, od.n)
    vec[1] = one(N)
    vec[2] = vec[1]
    return (vec, 1, 2)
end
function Base.next(od::OctDirections{N}, state::Tuple) where N
    if length(state) == 1
        # 1D case
        if state == 1
            return (ones(N, 1), 2)
        else
            @assert state == 2
            return (-ones(N, 1), 0)
        end
    end
    vec = state[1]
    i = state[2]
    j = state[3]
    if vec[j] > 0
        # change sign at j
        vec[j] = -one(N)
    elseif vec[i] > 0
        # change sign at i and j
        vec[i] = -one(N)
        vec[j] = one(N)
    elseif j < od.n
        # advance j, reset i
        vec[i] = one(N)
        vec[j] = zero(N)
        j = j + 1
        vec[j] = one(N)
    elseif i < od.n - 1
        # advance i
        vec[i] = zero(N)
        vec[j] = vec[i]
        i = i + 1
        j = i + 1
        vec[i] = one(N)
        vec[j] = vec[i]
    else
        vec = zeros(N, od.n)
        vec[1] = one(N)
        vec[2] = vec[1]
        return (vec, 1) # continue with box directions
    end
    return (copy(vec), (vec, i, j))
end
Base.next(od::OctDirections{N}, state::Int) where N =
    next(BoxDirections{N}(od.n), state)
Base.done(od::OctDirections, state) = state == 0
Base.length(od::OctDirections) = 2 * od.n^2

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
Base.next(bdd::BoxDiagDirections{N}, state::Int) where N =
    next(BoxDirections{N}(bdd.n), state)
Base.done(bdd::BoxDiagDirections, state) = state == 0
Base.length(bdd::BoxDiagDirections) = bdd.n == 1 ? 2 : 2^bdd.n + 2 * bdd.n
