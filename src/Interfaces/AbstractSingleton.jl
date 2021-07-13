import Base.∈

export AbstractSingleton,
       element,
       an_element,
       linear_map

"""
    AbstractSingleton{N} <: AbstractHyperrectangle{N}

Abstract type for sets with a single value.

### Notes

Every concrete `AbstractSingleton` must define the following functions:
- `element(::AbstractSingleton{N})` -- return the single element
- `element(::AbstractSingleton{N}, i::Int)` -- return the single element's
    entry in the `i`-th dimension

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractSingleton)
2-element Vector{Any}:
 Singleton
 ZeroSet
```
"""
abstract type AbstractSingleton{N} <: AbstractHyperrectangle{N} end

isconvextype(::Type{<:AbstractSingleton}) = true

# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(S::AbstractSingleton{N}, i::Int) where {N}

Return the box radius of a set with a single value in a given dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

Zero.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}, i::Int) where {N}
    return zero(N)
end


"""
    radius_hyperrectangle(S::AbstractSingleton{N}) where {N}

Return the box radius of a set with a single value in every dimension.

### Input

- `S` -- set with a single value

### Output

The zero vector.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}) where {N}
    return zeros(N, dim(S))
end


"""
    high(S::AbstractSingleton)

Return the higher coordinates of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A vector with the higher coordinates of the set with a single value.
"""
function high(S::AbstractSingleton)
    return element(S)
end

"""
    high(S::AbstractSingleton, i::Int)

Return the higher coordinate of a set with a single value in the given
dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

The higher coordinate of the set with a single value in the given dimension.
"""
function high(S::AbstractSingleton, i::Int)
    return element(S)[i]
end

"""
    low(S::AbstractSingleton)

Return the lower coordinates of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A vector with the lower coordinates of the set with a single value.
"""
function low(S::AbstractSingleton)
    return element(S)
end

"""
    low(S::AbstractSingleton, i::Int)

Return the lower coordinate of a set with a single value in the given
dimension.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

The lower coordinate of the set with a single value in the given dimension.
"""
function low(S::AbstractSingleton, i::Int)
    return element(S)[i]
end


# --- AbstractZonotope interface functions ---


"""
   genmat(S::AbstractSingleton{N}) where {N}

Return the (empty) generator matrix of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A matrix with no columns representing the generators of `S`.
"""
function genmat(S::AbstractSingleton{N}) where {N}
    return Matrix{N}(undef, dim(S), 0)
end

"""
    generators(S::AbstractSingleton{N}) where {N}

Return an (empty) iterator over the generators of a set with a single value.

### Input

- `S` -- set with a single value

### Output

An empty iterator.
"""
function generators(S::AbstractSingleton{N}) where {N}
    return EmptyIterator{Vector{N}}()
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
    center(S::AbstractSingleton)

Return the center of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element of the set.
"""
function center(S::AbstractSingleton)
    return element(S)
end


# --- AbstractPolytope interface functions ---


"""
    vertices(S::AbstractSingleton{N}) where {N}

Construct an iterator over the vertices of a set with a single value.

### Input

- `S` -- set with a single value

### Output

An iterator with a single value.
"""
function vertices(S::AbstractSingleton{N}) where {N}
    return SingletonIterator(element(S))
end

"""
    vertices_list(S::AbstractSingleton)

Return the list of vertices of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A list containing only a single vertex.
"""
function vertices_list(S::AbstractSingleton)
    return [element(S)]
end

"""
    linear_map(M::AbstractMatrix, S::AbstractSingleton)

Concrete linear map of an abstract singleton.

### Input

- `M` -- matrix
- `S` -- abstract singleton

### Output

The abstract singleton of the same type of ``S`` obtained by applying the
linear map to the element in ``S``.
"""
function linear_map(M::AbstractMatrix, S::AbstractSingleton)
    @assert dim(S) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(S))"

    T = typeof(S)
    return T(M * element(S))
end

# --- LazySet interface functions ---


"""
    σ(d::AbstractVector, S::AbstractSingleton)

Return the support vector of a set with a single value.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

The support vector, which is the set's vector itself, irrespective of the given
direction.
"""
function σ(d::AbstractVector, S::AbstractSingleton)
    return element(S)
end

"""
    ρ(d::AbstractVector, S::AbstractSingleton)

Evaluate the support function of a set with a single value in a given direction.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, S::AbstractSingleton)
    return dot(d, element(S))
end

"""
    ∈(x::AbstractVector, S::AbstractSingleton)

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
function ∈(x::AbstractVector, S::AbstractSingleton)
    return x == element(S)
end

# this operation is forbidden, but it is a common error
function ∈(S::AbstractSingleton, X::LazySet)
    throw(ArgumentError("cannot make a point-in-set check if the left-hand side is " *
          "a set; either check for set inclusion, as in `S ⊆ X`, or check for " *
          "membership, as in `element(S) ∈ X` (the results are equivalent but " *
          "the implementations may differ)"))
end

function chebyshev_center(S::AbstractSingleton{N}; compute_radius::Bool=false) where {N}
    if compute_radius
        return element(S), zero(N)
    end
    return element(S)
end

"""
    plot_recipe(S::AbstractSingleton{N}, [ε]=zero(N)) where {N}

Convert a singleton to a pair `(x, y)` of points for plotting.

### Input

- `S` -- singleton
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of one point that can be plotted.
"""
function plot_recipe(S::AbstractSingleton{N}, ε=zero(N)) where {N}
    @assert dim(S) <= 3 "cannot plot a $(dim(S))-dimensional $(typeof(S))"

    if dim(S) == 1
        return [element(S)[1]], [zero(N)]
    else
        return [element(S)[1]], [element(S)[2]]
    end
end
