import Base.rand

export Singleton

"""
    Singleton{N<:Real, VN<:AbstractVector{N}} <: AbstractSingleton{N}

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton{N<:Real, VN<:AbstractVector{N}} <: AbstractSingleton{N}
    element::VN
end

isoperationtype(::Type{<:Singleton}) = false

# --- AbstractSingleton interface functions ---


"""
    element(S::Singleton{N}) where {N<:Real}

Return the element of a singleton.

### Input

- `S` -- singleton

### Output

The element of the singleton.
"""
function element(S::Singleton{N}) where {N<:Real}
    return S.element
end

"""
    element(S::Singleton{N}, i::Int)::N where {N<:Real}

Return the i-th entry of the element of a singleton.

### Input

- `S` -- singleton
- `i` -- dimension

### Output

The i-th entry of the element of the singleton.
"""
function element(S::Singleton{N}, i::Int)::N where {N<:Real}
    return S.element[i]
end


# --- LazySet interface functions ---


"""
    rand(::Type{Singleton}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Singleton{N}

Create a random singleton.

### Input

- `Singleton` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random singleton.

### Algorithm

The element is a normally distributed vector with entries of mean 0 and standard
deviation 1.
"""
function rand(::Type{Singleton};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Singleton{N}
    rng = reseed(rng, seed)
    element = randn(rng, N, dim)
    return Singleton(element)
end

"""
    translate(S::Singleton{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) a singleton by a given vector.

### Input

- `S` -- singleton
- `v` -- translation vector

### Output

A translated singleton.

### Algorithm

We add the vector to the point in the singleton.
"""
function translate(S::Singleton{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(S) "cannot translate a $(dim(S))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Singleton(element(S) + v)
end
