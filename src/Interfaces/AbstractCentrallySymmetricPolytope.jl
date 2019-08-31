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
  - `center(::AbstractCentrallySymmetricPolytope{N})::Vector{N}` -- return the
     center point
- from `AbstractPolytope`:
  - `vertices_list(::AbstractCentrallySymmetricPolytope{N})::Vector{Vector{N}}`
     -- return a list of all vertices

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetricPolytope)
2-element Array{Any,1}:
 AbstractZonotope
 Ball1
```
"""
abstract type AbstractCentrallySymmetricPolytope{N<:Real} <: AbstractPolytope{N}
end


# --- common AbstractCentrallySymmetric functions (copy-pasted) ---


"""
    dim(P::AbstractCentrallySymmetricPolytope)::Int

Return the ambient dimension of a centrally symmetric, polytopic set.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

The ambient dimension of the polytopic set.
"""
@inline function dim(P::AbstractCentrallySymmetricPolytope)::Int
    return length(center(P))
end


"""
    an_element(P::AbstractCentrallySymmetricPolytope{N})::Vector{N}
        where {N<:Real}

Return some element of a centrally symmetric polytope.

### Input

- `P` -- centrally symmetric polytope

### Output

The center of the centrally symmetric polytope.
"""
function an_element(P::AbstractCentrallySymmetricPolytope{N}
                   )::Vector{N} where {N<:Real}
    return center(P)
end

"""
    isempty(P::AbstractCentrallySymmetricPolytope)::Bool

Return if a centrally symmetric, polytopic set is empty or not.

### Input

- `P` -- centrally symmetric, polytopic set

### Output

`false`.
"""
function isempty(P::AbstractCentrallySymmetricPolytope)::Bool
    return false
end
