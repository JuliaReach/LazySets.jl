export AbstractPolyhedron

"""
    AbstractPolyhedron{N<:Real} <: LazySet{N}

Abstract type for polyhedral sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces.

### Notes

Every concrete `AbstractPolyhedron` must define the following functions:
- `constraints_list(::AbstractPolyhedron{N})::Vector{LinearConstraint{N}}` --
    return a list of all facet constraints

```jldoctest
julia> subtypes(AbstractPolyhedron)
5-element Array{Any,1}:
 AbstractPolytope
 HPolyhedron
 HalfSpace
 Hyperplane
 Line
```
"""
abstract type AbstractPolyhedron{N<:Real} <: LazySet{N} end


# --- common AbstractPolyhedron functions ---


# To account for the compilation order, functions are defined in the file
# AbstractPolyhedron_functions.jl
