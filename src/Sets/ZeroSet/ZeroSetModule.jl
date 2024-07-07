module ZeroSetModule

using Reexport

using ..LazySets: AbstractSingleton, Singleton
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: dim, isoperationtype, rand, rectify, reflect, ∈,
                        linear_map, scale, scale!, ρ, translate
@reexport import ..LazySets: element
@reexport using ..API

export ZeroSet

"""
    ZeroSet{N} <: AbstractSingleton{N}

Type that represents the zero set, i.e., the set that only contains the origin.

### Fields

- `dim` -- the ambient dimension of the set
"""
struct ZeroSet{N} <: AbstractSingleton{N}
    dim::Int
end

isoperationtype(::Type{<:ZeroSet}) = false

# default constructor of type Float64
ZeroSet(dim::Int) = ZeroSet{Float64}(dim)

"""
    element(Z::ZeroSet{N}) where {N}

Return the element of a zero set.

### Input

- `Z` -- zero set

### Output

The element of the zero set, i.e., a zero vector.
"""
function element(Z::ZeroSet{N}) where {N}
    return zeros(N, Z.dim)
end

"""
    element(Z::ZeroSet{N}, ::Int) where {N}

Return the i-th entry of the element of a zero set.

### Input

- `Z` -- zero set
- `i` -- dimension

### Output

The i-th entry of the element of the zero set, i.e., 0.
"""
function element(Z::ZeroSet{N}, ::Int) where {N}
    return zero(N)
end

"""
    dim(Z::ZeroSet)

Return the ambient dimension of a zero set.

### Input

- `Z` -- zero set

### Output

The ambient dimension of the zero set.
"""
function dim(Z::ZeroSet)
    return Z.dim
end

"""
    ρ(d::AbstractVector, Z::ZeroSet)

Evaluate the support function of a zero set in a given direction.

### Input

- `d` -- direction
- `Z` -- zero set

### Output

`0`.
"""
function ρ(d::AbstractVector, Z::ZeroSet)
    @assert length(d) == dim(Z) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(Z))-dimensional set"
    N = promote_type(eltype(d), eltype(Z))
    return zero(N)
end

"""
    ∈(x::AbstractVector, Z::ZeroSet)

Check whether a given point is contained in a zero set.

### Input

- `x` -- point/vector
- `Z` -- zero set

### Output

`true` iff ``x ∈ Z``.

### Examples

```jldoctest
julia> Z = ZeroSet(2);

julia> [1.0, 0.0] ∈ Z
false
julia> [0.0, 0.0] ∈ Z
true
```
"""
function ∈(x::AbstractVector, Z::ZeroSet)
    @assert length(x) == dim(Z) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(Z))-dimensional set"
    return iszero(x)
end

"""
    rand(::Type{ZeroSet}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a zero set (note that there is nothing to randomize).

### Input

- `ZeroSet` -- type for dispatch
- `N`       -- (optional, default: `Float64`) numeric type
- `dim`     -- (optional, default: 2) dimension
- `rng`     -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`    -- (optional, default: `nothing`) seed for reseeding

### Output

The (only) zero set of the given numeric type and dimension.
"""
function rand(::Type{ZeroSet};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return ZeroSet{N}(dim)
end

"""
    linear_map(M::AbstractMatrix, Z::ZeroSet)

Concrete linear map of a zero set.

### Input

- `M` -- matrix
- `Z` -- zero set

### Output

The zero set whose dimension matches the output dimension of the given matrix.
"""
function linear_map(M::AbstractMatrix, Z::ZeroSet)
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"

    N = promote_type(eltype(M), eltype(Z))
    return ZeroSet{N}(size(M, 1))
end

"""
    translate(Z::ZeroSet, v::AbstractVector)

Translate (i.e., shift) a zero set by a given vector.

### Input

- `Z` -- zero set
- `v` -- translation vector

### Output

A singleton containing the vector `v`.
"""
function translate(Z::ZeroSet, v::AbstractVector)
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Singleton(v)
end

"""
    rectify(Z::ZeroSet)

Concrete rectification of a zero set.

### Input

- `Z` -- zero set

### Output

The same set.
"""
function rectify(Z::ZeroSet)
    return Z
end

"""
    reflect(Z::ZeroSet)

Concrete reflection of a zero set `Z`, resulting in the reflected set `-Z`.

### Input

- `Z` -- zero set

### Output

The same zero set.
"""
function reflect(Z::ZeroSet)
    return Z
end

function scale(::Real, Z::ZeroSet)
    return Z
end

function scale!(::Real, Z::ZeroSet)
    return Z
end

end  # module
