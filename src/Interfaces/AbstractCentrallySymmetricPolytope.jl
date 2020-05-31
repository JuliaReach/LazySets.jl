import Base.isempty

export AbstractCentrallySymmetricPolytope,
       center,
       an_element,
       vertices_list,
       singleton_list

"""
    AbstractCentrallySymmetricPolytope{N<:Real} <: AbstractPolytope{N}

Abstract type for centrally symmetric, polytopic sets.
It combines the `AbstractCentrallySymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractCentrallySymmetricPolytope` must define the following
functions:
- from `AbstractCentrallySymmetric`:
  - `center(::AbstractCentrallySymmetricPolytope{N})` -- return the
     center point
- from `AbstractPolytope`:
  - `vertices_list(::AbstractCentrallySymmetricPolytope{N})`
     -- return a list of all vertices

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetricPolytope)
2-element Array{Any,1}:
 AbstractZonotope
 Ball1
```
"""
abstract type AbstractCentrallySymmetricPolytope{N<:Real} <: AbstractPolytope{N} end

isconvextype(::Type{<:AbstractCentrallySymmetricPolytope}) = true


# --- common AbstractCentrallySymmetric functions (copy-pasted) ---


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
    an_element(P::AbstractCentrallySymmetricPolytope{N}) where {N<:Real}

Return some element of a centrally symmetric polytope.

### Input

- `P` -- centrally symmetric polytope

### Output

The center of the centrally symmetric polytope.
"""
function an_element(P::AbstractCentrallySymmetricPolytope{N}) where {N<:Real}
    return center(P)
end

"""
    isempty(P::AbstractCentrallySymmetricPolytope)

Return if a centrally symmetric, polytopic set is empty or not.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

`false`.
"""
function isempty(::AbstractCentrallySymmetricPolytope)
    return false
end

"""
    isuniversal(S::AbstractCentrallySymmetricPolytope{N}, [witness]::Bool=false
               ) where {N<:Real}

Check whether a centrally symmetric polytope is universal.

### Input

- `S`       -- centrally symmetric polytope
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``

### Algorithm

A witness is obtained by computing the support vector in direction
`d = [1, 0, …, 0]` and adding `d` on top.
"""
function isuniversal(S::AbstractCentrallySymmetricPolytope{N},
                     witness::Bool=false) where {N<:Real}
    if witness
        d = SingleEntryVector{N}(1, dim(S))
        w = σ(d, S) + d
        return (false, w)
    else
        return false
    end
end

"""
    center(H::AbstractCentrallySymmetricPolytope{N}, i::Int) where {N<:Real}

Return the center along a given dimension of a centrally symmetric polytope.

### Input

- `S` -- centrally symmetric polytope
- `i` -- dimension of interest

### Output

The center along a given dimension of the centrally symmetric polytope.
"""
@inline function center(S::AbstractCentrallySymmetricPolytope{N},
                        i::Int) where {N<:Real}
    return center(S)[i]
end
