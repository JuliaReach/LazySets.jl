export AbstractCentrallySymmetric

"""
    AbstractCentrallySymmetric{N} <: ConvexSet{N}

Abstract type for centrally symmetric compact convex sets.

### Notes

Every concrete `AbstractCentrallySymmetric` must define the following function:

- `center(::AbstractCentrallySymmetric)` -- return the center point

The subtypes of `AbstractCentrallySymmetric`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetric)
2-element Vector{Any}:
 AbstractBallp
 Ellipsoid
```
"""
abstract type AbstractCentrallySymmetric{N} <: ConvexSet{N} end

# a set with a unique center must be bounded
function isboundedtype(::Type{<:AbstractCentrallySymmetric})
    return true
end

# a set with a unique center must be bounded
function isbounded(::AbstractCentrallySymmetric)
    return true
end

# To account for the compilation order and sharing with
# AbstractCentrallySymmetricPolytope, other functions are defined in the file
# AbstractCentrallySymmetric_functions.jl
# The functions above do not need to be shared.
