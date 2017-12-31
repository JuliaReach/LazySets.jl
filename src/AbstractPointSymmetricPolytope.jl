export AbstractPointSymmetricPolytope,
       center,
       an_element,
       vertices_list,
       singleton_list

"""
    AbstractPointSymmetricPolytope{N<:Real} <: AbstractPolytope{N}

Abstract type for point symmetric, polytopic sets.
It combines the `AbstractPointSymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractPointSymmetricPolytope` must define the following
functions:
- from `AbstractPointSymmetric`:
  - `center(::AbstractPointSymmetricPolytope{N})::Vector{N}` -- return the
     center point
- from `AbstractPolytope`:
  - `vertices_list(::AbstractPointSymmetricPolytope{N})::Vector{Vector{N}}`
     -- return a list of all vertices

```jldoctest
julia> subtypes(AbstractPointSymmetricPolytope)
3-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractHyperrectangle
 LazySets.Ball1
 LazySets.Zonotope
```
"""
abstract type AbstractPointSymmetricPolytope{N<:Real} <: AbstractPolytope{N} end


# --- common AbstractPointSymmetric functions (copy-pasted) ---


"""
    dim(P::AbstractPointSymmetricPolytope)::Int

Return the ambient dimension of a point symmetric set.

### Input

- `P` -- set

### Output

The ambient dimension of the set.
"""
@inline function dim(P::AbstractPointSymmetricPolytope)::Int
    return length(center(P))
end


"""
    an_element(P::AbstractPointSymmetricPolytope{N})::Vector{N} where {N<:Real}

Return some element of a point symmetric polytope.

### Input

- `P` -- point symmetric polytope

### Output

The center of the point symmetric polytope.
"""
function an_element(P::AbstractPointSymmetricPolytope{N}
                   )::Vector{N} where {N<:Real}
    return center(P)
end
