import Base.∈

export AbstractSingleton,
       element,
       an_element,
       linear_map

"""
    AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N}

Abstract type for sets with a single value.

### Notes

Every concrete `AbstractSingleton` must define the following functions:
- `element(::AbstractSingleton{N})::Vector{N}` -- return the single element
- `element(::AbstractSingleton{N}, i::Int)::N` -- return the single element's
    entry in the `i`-th dimension

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractSingleton)
2-element Array{Any,1}:
 Singleton
 ZeroSet
```
"""
abstract type AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N} end

isconvextype(::Type{<:AbstractSingleton}) = true

# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}

Return the box radius of a set with a single value in a given dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

Zero.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}, i::Int
                              )::N where {N<:Real}
    return zero(N)
end


"""
    radius_hyperrectangle(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the box radius of a set with a single value in every dimension.

### Input

- `S` -- set with a single value

### Output

The zero vector.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}
                              )::Vector{N} where {N<:Real}
    return zeros(N, dim(S))
end


"""
    high(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the higher coordinates of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A vector with the higher coordinates of the set with a single value.
"""
function high(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end

"""
    high(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}

Return the higher coordinate of a set with a single value in the given
dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

The higher coordinate of the set with a single value in the given dimension.
"""
function high(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}
    return element(S)[i]
end

"""
    low(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the lower coordinates of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A vector with the lower coordinates of the set with a single value.
"""
function low(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end

"""
    low(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}

Return the lower coordinate of a set with a single value in the given
dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

The lower coordinate of the set with a single value in the given dimension.
"""
function low(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}
    return element(S)[i]
end


# --- AbstractZonotope interface functions ---


"""
   genmat(S::AbstractSingleton)

Return the (empty) generator matrix of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A matrix with no columns representing the generators of `S`.
"""
function genmat(S::AbstractSingleton{N}) where {N<:Real}
    return Matrix{N}(undef, dim(S), 0)
end

# iterator that wraps the generator matrix
struct EmptyGeneratorIterator{N<:Real}
end

Base.length(::EmptyGeneratorIterator) = 0

Base.eltype(::Type{EmptyGeneratorIterator{N}}) where {N} = Vector{N}

function Base.iterate(::EmptyGeneratorIterator, state=nothing)
    return nothing
end

"""
    generators(S::AbstractSingleton)

Return an (empty) iterator over the generators of a set with a single value.

### Input

- `S` -- set with a single value

### Output

An empty iterator.
"""
function generators(S::AbstractSingleton{N}) where {N<:Real}
    return EmptyGeneratorIterator{N}()
end

"""
    ngens(S::AbstractSingleton)

Return the number of generators of a set with a single value.

### Input

- `H` -- set with a single value

### Output

The number of generators.

### Algorithm

A set with a single value has no generators, so the result is ``0``.
"""
function ngens(S::AbstractSingleton)
    return 0
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the center of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element of the set.
"""
function center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(S::AbstractSingleton{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A list containing only a single vertex.
"""
function vertices_list(S::AbstractSingleton{N}
                      )::Vector{Vector{N}} where {N<:Real}
    return [element(S)]
end

"""
    linear_map(M::AbstractMatrix{N}, S::AbstractSingleton{N}) where {N<:Real}

Concrete linear map of an abstract singleton.

### Input

- `M` -- matrix
- `S` -- abstract singleton

### Output

The abstract singleton of the same type of ``S`` obtained by applying the
linear map to the element in ``S``.
"""
function linear_map(M::AbstractMatrix{N},
                    S::AbstractSingleton{N}) where {N<:Real}
    @assert dim(S) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(S))"

    T = typeof(S)
    return T(M * element(S))
end

# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, S::AbstractSingleton{N}) where {N<:Real}

Return the support vector of a set with a single value.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

The support vector, which is the set's vector itself, irrespective of the given
direction.
"""
function σ(d::AbstractVector{N}, S::AbstractSingleton{N}) where {N<:Real}
    return element(S)
end

"""
    ρ(d::AbstractVector{N}, S::AbstractSingleton{N}) where {N<:Real}

Evaluate the support function of a set with a single value in a given direction.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector{N}, S::AbstractSingleton{N}) where {N<:Real}
    return dot(d, element(S))
end

"""
    ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}

Check whether a given point is contained in a set with a single value.

### Input

- `x` -- point/vector
- `S` -- set with a single value

### Output

`true` iff ``x ∈ S``.

### Notes

This implementation performs an exact comparison, which may be insufficient with
floating point computations.
"""
function ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}
    return x == element(S)
end

# this operation is forbidden, but it is a common error
function ∈(S::AbstractSingleton{N}, X::LazySet{N})::Bool where {N<:Real}
    error("cannot make a point-in-set check if the left-hand side is " *
          "a set; either check for set inclusion, as in `S ⊆ X`, or check for " *
          "membership, as in `element(S) ∈ X` (the results are equivalent but " *
          "the implementations may differ)")
end

"""
    plot_recipe(S::AbstractSingleton{N}, [ε]::N=zero(N)) where {N<:Real}

Convert a singleton to a pair `(x, y)` of points for plotting.

### Input

- `S` -- singleton
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of one point that can be plotted.
"""
function plot_recipe(S::AbstractSingleton{N}, ε::N=zero(N)) where {N<:Real}
    @assert dim(S) <= 3 "cannot plot a $(dim(S))-dimensional $(typeof(S))"

    if dim(S) == 1
        return [element(S)[1]], [zero(N)]
    else
        return [element(S)[1]], [element(S)[2]]
    end
end
