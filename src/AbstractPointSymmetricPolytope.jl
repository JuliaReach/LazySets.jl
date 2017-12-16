export AbstractPointSymmetricPolytope,
       center,
       vertices_list,
       singleton_list

"""
    AbstractPointSymmetricPolytope{N<:Real} <: LazySet

Abstract type for point symmetric, polytopic sets.
It combines the `AbstractPointSymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractPointSymmetricPolytope` must define the following
functions:
- from `AbstractPointSymmetric`:

  - `center(::AbstractPointSymmetric{N})::Vector{N}` -- return the center point
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
abstract type AbstractPointSymmetricPolytope{N<:Real} <: LazySet end


# --- common AbstractPointSymmetric functions (copy-pasted) ---


"""
    dim(S::AbstractPointSymmetricPolytope)::Int

Return the ambient dimension of a point symmetric set.

### Input

- `S` -- set

### Output

The ambient dimension of the set.
"""
@inline function dim(S::AbstractPointSymmetricPolytope)::Int
    return length(center(S))
end


# --- common AbstractPolytope functions (copy-pasted) ---


"""
    singleton_list(P::AbstractPointSymmetricPolytope{N}
                  )::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polytopic as a list of singletons.

### Input

- `P` -- a polytopic set

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::AbstractPointSymmetricPolytope{N}
                       )::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in P.vertices_list]
end
