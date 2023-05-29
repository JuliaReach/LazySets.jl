export AbstractPolyhedron

"""
    AbstractPolyhedron{N} <: ConvexSet{N}

Abstract type for compact convex polyhedral sets.

### Notes

Every concrete `AbstractPolyhedron` must define the following functions:

- `constraints_list(::AbstractPolyhedron{N})` -- return a list of all facet
    constraints

Polyhedra are defined as the intersection of a finite number of closed
half-spaces.
As such, polyhedra are closed and convex but not necessarily bounded.
Bounded polyhedra are called *polytopes* (see [`AbstractPolytope`](@ref)).

The subtypes of `AbstractPolyhedron` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolyhedron)
8-element Vector{Any}:
 AbstractPolytope
 HPolyhedron
 HalfSpace
 Hyperplane
 Line
 Line2D
 Star
 Universe
```
"""
abstract type AbstractPolyhedron{N} <: ConvexSet{N} end

# supertype for all linear map algorithms
abstract type AbstractLinearMapAlgorithm end

# To account for the compilation order, functions are defined in the file
# AbstractPolyhedron_functions.jl
