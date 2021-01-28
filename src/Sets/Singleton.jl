import Base.rand

export Singleton

"""
    Singleton{N, VN<:AbstractVector{N}} <: AbstractSingleton{N}

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton{N, VN<:AbstractVector{N}} <: AbstractSingleton{N}
    element::VN
end

isoperationtype(::Type{<:Singleton}) = false
isconvextype(::Type{<:Singleton}) = true

# --- AbstractSingleton interface functions ---


"""
    element(S::Singleton)

Return the element of a singleton.

### Input

- `S` -- singleton

### Output

The element of the singleton.
"""
function element(S::Singleton)
    return S.element
end

"""
    element(S::Singleton, i::Int)

Return the i-th entry of the element of a singleton.

### Input

- `S` -- singleton
- `i` -- dimension

### Output

The i-th entry of the element of the singleton.
"""
function element(S::Singleton, i::Int)
    return S.element[i]
end


# --- LazySet interface functions ---


"""
    rand(::Type{Singleton}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    element = randn(rng, N, dim)
    return Singleton(element)
end

"""
    translate(S::Singleton, v::AbstractVector)

Translate (i.e., shift) a singleton by a given vector.

### Input

- `S` -- singleton
- `v` -- translation vector

### Output

A translated singleton.

### Algorithm

We add the vector to the point in the singleton.
"""
function translate(S::Singleton, v::AbstractVector)
    @assert length(v) == dim(S) "cannot translate a $(dim(S))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Singleton(element(S) + v)
end

"""
    rectify(S::Singleton)

Concrete rectification of a singleton.

### Input

- `S` -- singleton

### Output

The `Singleton` that corresponds to the rectification of `S`.
"""
function rectify(S::Singleton)
    return Singleton(rectify(element(S)))
end

"""
    project(S::Singleton, block::AbstractVector{Int})

Concrete projection of a singleton.

### Input

- `S`     -- singleton
- `block` -- block structure, a vector with the dimensions of interest

### Output

A set representing the projection of the singleton `S` on the dimensions
specified by `block`.
"""
function project(S::Singleton, block::AbstractVector{Int})
    return Singleton(element(S)[block])
end
