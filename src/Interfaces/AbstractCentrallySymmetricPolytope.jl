export AbstractCentrallySymmetricPolytope

"""
    AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N}

Abstract type for centrally symmetric, polytopic sets.
It combines the `AbstractCentrallySymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractCentrallySymmetricPolytope` must define the following
methods:

- from `AbstractCentrallySymmetric`:
  - `center(::AbstractCentrallySymmetricPolytope)` -- return the center point
  - `center(::AbstractCentrallySymmetricPolytope, i::Int)` -- return the center
                                                              point at index `i`
- from `AbstractPolytope`:
  - `vertices_list(::AbstractCentrallySymmetricPolytope)` -- return a list of
                                                             all vertices

The subtypes of `AbstractCentrallySymmetricPolytope` (including abstract
interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetricPolytope)
2-element Vector{Any}:
 AbstractZonotope
 Ball1
```
"""
abstract type AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N} end

# common AbstractCentrallySymmetric functions
# copy-pasted because Julia does not have multiple inheritance

"""
    dim(P::AbstractCentrallySymmetricPolytope)

Return the ambient dimension of a centrally symmetric, polytopic set.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

The ambient dimension of the polytopic set.
"""
@inline function dim(P::AbstractCentrallySymmetricPolytope)
    return length(center(P))
end

"""
    an_element(P::AbstractCentrallySymmetricPolytope)

Return some element of a centrally symmetric, polytopic set.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

The center of the centrally symmetric, polytopic set.
"""
function an_element(P::AbstractCentrallySymmetricPolytope)
    return center(P)
end

"""
    isempty(P::AbstractCentrallySymmetricPolytope)

Check whether a centrally symmetric, polytopic set is empty.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

`false`.
"""
function isempty(::AbstractCentrallySymmetricPolytope)
    return false
end

"""
    isuniversal(S::AbstractCentrallySymmetricPolytope{N},
                [witness]::Bool=false) where {N}

Check whether a centrally symmetric, polytopic set is universal.

### Input

- `S`       -- centrally symmetric, polytopic set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``

### Algorithm

Centrally symmetric, polytopic sets are bounded.
A witness is obtained by computing the support vector in direction
`d = [1, 0, …, 0]` and adding `d` on top.
"""
function isuniversal(S::AbstractCentrallySymmetricPolytope{N},
                     witness::Bool=false) where {N}
    if witness
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end

"""
    center(S::AbstractCentrallySymmetricPolytope, i::Int)

Return the center of a centrally symmetric, polytopic set along a given
dimension.

### Input

- `S` -- centrally symmetric, polytopic set
- `i` -- dimension of interest

### Output

The center along the given dimension.
"""
@inline function center(S::AbstractCentrallySymmetricPolytope, i::Int)
    return center(S)[i]
end

"""
    extrema(S::AbstractCentrallySymmetricPolytope)

Return two vectors with the lowest and highest coordinate of a centrally
symmetric, polytopic set.

### Input

- `S` -- centrally symmetric, polytopic set

### Output

Two vectors with the lowest and highest coordinates of `S`.

### Notes

The result is equivalent to `(low(S), high(S))`.

### Algorithm

We compute `high(S)` and then compute the lowest coordinates with the help of
`center(S)` (which is assumed to be cheaper to obtain).
"""
function extrema(S::AbstractCentrallySymmetricPolytope)
    # h = c + r
    h = high(S)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 .* center(S) .- h
    return (l, h)
end

"""
    extrema(S::AbstractCentrallySymmetricPolytope, i::Int)

Return the lower and higher coordinate of a centrally symmetric, polytopic set
in a given dimension.

### Input

- `S` -- centrally symmetric, polytopic set
- `i` -- dimension of interest

### Output

The lower and higher coordinate of the centrally symmetric, polytopic set in the
given dimension.

### Notes

The result is equivalent to `(low(S, i), high(S, i))`.

### Algorithm

We compute `high(S, i)` and then compute the lowest coordinates with the help of
`center(S, i)` (which is assumed to be cheaper to obtain).
"""
function extrema(S::AbstractCentrallySymmetricPolytope, i::Int)
    # h = c + r
    h = high(S, i)
    # l = c - r = -c - r + 2 * c = 2 * c - h
    l = 2 * center(S, i) - h
    return (l, h)
end
