export AbstractSingleton,
       element

"""
    AbstractSingleton{N} <: AbstractHyperrectangle{N}

Abstract type for sets with a single value.

### Notes

See [`Singleton`](@ref) for a standard implementation of this interface.

Every concrete `AbstractSingleton` must define the following function:

- `element(::AbstractSingleton)` -- return the single element

Among other functions, the following function is then automatically defined:

- `element(::AbstractSingleton, i::Int)` -- return the single element at index
                                            `i`

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractSingleton)
2-element Vector{Any}:
 Singleton
 ZeroSet
```
"""
abstract type AbstractSingleton{N} <: AbstractHyperrectangle{N} end

"""
    element(S::AbstractSingleton)

Return the element of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The unique element of `S`.
"""
function element(::AbstractSingleton) end

"""
    element(S::AbstractSingleton, i::Int)

Return the i-th entry of the element of a set with a single value.

### Input

- `S` -- set with a single value
- `i` -- dimension of interest

### Output

The i-th entry of the element.
"""
function element(S::AbstractSingleton, i::Int)
    return element(S)[i]
end

@validate function radius_hyperrectangle(S::AbstractSingleton, i::Int)
    N = eltype(S)
    return zero(N)
end

function radius_hyperrectangle(S::AbstractSingleton)
    N = eltype(S)
    return zeros(N, dim(S))
end

function high(S::AbstractSingleton)
    return element(S)
end

@validate function high(S::AbstractSingleton, i::Int)
    return element(S)[i]
end

function low(S::AbstractSingleton)
    return element(S)
end

@validate function low(S::AbstractSingleton, i::Int)
    return element(S)[i]
end

function genmat(S::AbstractSingleton)
    N = eltype(S)
    return Matrix{N}(undef, dim(S), 0)
end

function generators(S::AbstractSingleton)
    N = eltype(S)
    return EmptyIterator{Vector{N}}()
end

function ngens(::AbstractSingleton)
    return 0
end

function center(S::AbstractSingleton)
    return element(S)
end

@validate function center(S::AbstractSingleton, i::Int)
    return element(S, i)
end

function vertices(S::AbstractSingleton)
    return SingletonIterator(element(S))
end

function vertices_list(S::AbstractSingleton)
    return [element(S)]
end

"""
# Extended help

    σ(d::AbstractVector, S::AbstractSingleton)

### Algorithm

The support vector is the set's vector itself, irrespective of the given
direction.
"""
@validate function σ(d::AbstractVector, S::AbstractSingleton)
    return element(S)
end

@validate function ρ(d::AbstractVector, S::AbstractSingleton)
    return dot(d, element(S))
end

"""
# Extended help

    in(x::AbstractVector, S::AbstractSingleton)

### Notes

This implementation performs an approximate comparison to account for
imprecision in floating-point computations.
"""
@validate function in(x::AbstractVector, S::AbstractSingleton)
    return _isapprox(x, element(S))
end

# this operation is forbidden, but it is a common error
function in(S::AbstractSingleton, X::LazySet)
    throw(ArgumentError("cannot make a point-in-set check if the left-hand " *
                        "side is a set; either check for set inclusion, as in `S ⊆ X`, or " *
                        "check for membership, as in `element(S) ∈ X` (the results are " *
                        "equivalent, but the implementations may differ)"))
end

function chebyshev_center_radius(S::AbstractSingleton{N}) where {N}
    return element(S), zero(N)
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
    n = dim(S)
    if n == 1
        return [element(S)[1]], [zero(N)]
    elseif n == 2
        return [element(S)[1]], [element(S)[2]]
    elseif n == 3
        return [element(S)[1]], [element(S)[2]], [element(S)[3]]
    else
        throw(ArgumentError("cannot plot a $n-dimensional $(typeof(S))"))
    end
end

"""
# Extended help

    reflect(S::AbstractSingleton)

### Output

A `Singleton`.
"""
function reflect(S::AbstractSingleton)
    return Singleton(-element(S))
end

function constraints_list(S::AbstractSingleton; min_constraints::Bool=false)
    if min_constraints
        # fewest constraints (n+1) but more expensive to represent (`Vector`)
        return _constraints_list_singleton_Vector(element(S))
    else
        # more constraints (2n) but cheaper to represent (`SingleEntryVector`)
        return _constraints_list_hyperrectangle(S)
    end
end

# fewest constraints (n+1)
# Note: constraints are sorted CCW in 2D
function _constraints_list_singleton_Vector(e::AbstractVector{N}) where {N}
    n = length(e)
    constraints = Vector{HalfSpace{N,Vector{N}}}(undef, n + 1)
    @inbounds for i in 1:n
        # x_i >= e
        ai = zeros(N, n)
        ai[i] = -one(N)
        constraints[i] = HalfSpace(ai, -e[i])
    end
    # (∑_i x_i) <= ∑_i e_i
    @inbounds constraints[end] = HalfSpace(ones(N, n), sum(e))
    return constraints
end
