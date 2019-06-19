"""
    sign_cadlag(x::N)::N where {N<:Real}

This function works like the sign function but is ``1`` for input ``0``.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``-1`` otherwise.

### Notes

This is the sign function right-continuous at zero (see
[càdlàg function](https://en.wikipedia.org/wiki/C%C3%A0dl%C3%A0g)).
It can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> LazySets.sign_cadlag.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::N)::N where {N<:Real}
    return x < zero(x) ? -one(x) : one(x)
end

"""
    substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector

### Output

A fresh vector corresponding to `x` after `substitution` was applied.
"""
function substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(substitution, copy(x))
end

"""
    substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector (modified in this function)

### Output

The same (but see the Notes below) vector `x` but after `substitution` was
applied.

### Notes

The vector `x` is modified in-place if it has type `Vector` or `SparseVector`.
Otherwise, we first create a new `Vector` from it.
"""
function substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(Vector(x), substitution)
end

function substitute!(substitution::Dict{Int, T},
                     x::Union{Vector{T}, SparseVector{T}}) where {T}
    for (index, value) in substitution
        x[index] = value
    end
    return x
end

"""
    reseed(rng::AbstractRNG, seed::Union{Int, Nothing})::AbstractRNG

Reset the RNG seed if the seed argument is a number.

### Input

- `rng`  -- random number generator
- `seed` -- seed for reseeding

### Output

The input RNG if the seed is `nothing`, and a reseeded RNG otherwise.
"""
function reseed(rng::AbstractRNG, seed::Union{Int, Nothing})::AbstractRNG
    if seed != nothing
        return Random.seed!(rng, seed)
    end
    return rng
end

"""
    StrictlyIncreasingIndices

Iterator over the vectors of `m` strictly increasing indices from 1 to `n`.

### Fields

- `n` -- size of the index domain
- `m` -- number of indices to choose (resp. length of the vectors)

### Notes

The vectors are modified in-place.

The iterator ranges over ``\\binom{n}{m}`` (`n` choose `m`) possible vectors.

This implementation results in a lexicographic order with the last index growing
first.

### Examples

```jldoctest
julia> for v in LazySets.StrictlyIncreasingIndices(4, 2)
           println(v)
       end
[1, 2]
[1, 3]
[1, 4]
[2, 3]
[2, 4]
[3, 4]
```
"""
struct StrictlyIncreasingIndices
    n::Int
    m::Int

    function StrictlyIncreasingIndices(n::Int, m::Int)
        @assert n >= m > 0 "require n >= m > 0"
        new(n, m)
    end
end

Base.eltype(::Type{StrictlyIncreasingIndices}) = Vector{Int}
Base.length(sii::StrictlyIncreasingIndices) = binomial(sii.n, sii.m)

# initialization
function Base.iterate(sii::StrictlyIncreasingIndices)
    v = [1:sii.m;]
    return (v, v)
end

# normal iteration
function Base.iterate(sii::StrictlyIncreasingIndices, state::AbstractVector{Int})
    v = state
    i = sii.m
    diff = sii.n
    if i == diff
        return nothing
    end
    while v[i] == diff
        i -= 1
        diff -= 1
    end
    # update vector
    v[i] += 1
    for j in i+1:sii.m
        v[j] = v[j-1] + 1
    end
    # detect termination: first index has maximum value
    if i == 1 && v[1] == (sii.n - sii.m + 1)
        return (v, nothing)
    end
    return (v, v)
end

# termination
function Base.iterate(sii::StrictlyIncreasingIndices, state::Nothing)
    return nothing
end

"""
    subtypes(interface, concrete::Bool)

Return the concrete subtypes of a given interface.

### Input

- `interface` -- an abstract type, usually a set interface
- `concrete`  -- if `true`, seek further the inner abstract subtypes of the given
                 interface, otherwise return only the direct subtypes of `interface`

### Output

A list with the subtypes of the abstract type `interface`, sorted alphabetically.

### Examples

Consider the `AbstractPolytope` interface. If we include the abstract subtypes
of this interface,

```jldoctest subtypes
julia> using LazySets: subtypes

julia> subtypes(AbstractPolytope, false)
4-element Array{Any,1}:
 AbstractCentrallySymmetricPolytope
 AbstractPolygon
 HPolytope
 VPolytope
```

We can use this function to obtain the concrete subtypes of
`AbstractCentrallySymmetricPolytope` and `AbstractPolygon` (further until all
concrete types are obtained), using the `concrete` flag:

```jldoctest subtypes
julia> subtypes(AbstractPolytope, true)
14-element Array{Type,1}:
 Ball1
 BallInf
 HPolygon
 HPolygonOpt
 HPolytope
 Hyperrectangle
 Interval
 LineSegment
 Singleton
 SymmetricIntervalHull
 VPolygon
 VPolytope
 ZeroSet
 Zonotope
```
"""
function subtypes(interface, concrete::Bool)

    subtypes_to_test = subtypes(interface)

    # do not seek the concrete subtypes further
    if !concrete
        return sort(subtypes_to_test, by=string)
    end

    result = Vector{Type}()
    i = 0
    while i < length(subtypes_to_test)
        i += 1
        subtype = subtypes_to_test[i]
        new_subtypes = subtypes(subtype)
        if isempty(new_subtypes)
            # base type found
            push!(result, subtype)
        else
            # yet another interface layer
            append!(subtypes_to_test, new_subtypes)
        end
    end
    return sort(result, by=string)
end
