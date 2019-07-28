export AbstractPolyhedron

"""
    AbstractPolyhedron{N<:Real} <: LazySet{N}

Abstract type for compact convex polyhedral sets.

### Notes

Every concrete `AbstractPolyhedron` must define the following functions:
- `constraints_list(::AbstractPolyhedron{N})` -- return a list of all facet
    constraints

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolyhedron)
6-element Array{Any,1}:
 AbstractPolytope
 HPolyhedron
 HalfSpace
 Hyperplane
 Line
 Universe
```

Polyhedra are defined as the intersection of a finite number of closed
half-spaces.
As such, polyhedra are closed and convex but not necessarily bounded.
Bounded polyhedra are called *polytopes* (see [`AbstractPolytope`](@ref)).
"""
abstract type AbstractPolyhedron{N<:Real} <: LazySet{N} end


# --- common AbstractPolyhedron functions ---


# To account for the compilation order, functions are defined in the file
# AbstractPolyhedron_functions.jl
