"""
    IntersectionArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

This type assumes that the dimensions of all elements match.

The `EmptySet` is the absorbing element for `IntersectionArray`.

The intersection preserves convexity: if the set arguments are convex, then
their intersection is convex as well.

The convenience alias `∩` can be typed by `\\cap<tab>`.
"""
struct IntersectionArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    Intersection!(X, Y)

Convenience function to compute the lazy intersection and modify
`IntersectionArray`s in-place.
"""
function Intersection! end

∩(X::LazySet, Xs::LazySet...) = IntersectionArray(vcat(X, Xs...))
∩(X::LazySet) = X
∩(Xs::Vector{<:LazySet}) = IntersectionArray(Xs)

# constructor for an empty sum with optional size hint and numeric type
function IntersectionArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return IntersectionArray(arr)
end

# Universe is the neutral element for IntersectionArray
@neutral(IntersectionArray, Universe)

# EmptySet is the absorbing element for IntersectionArray
@absorbing(IntersectionArray, EmptySet)

# add functions connecting Intersection and IntersectionArray
@declare_array_version(Intersection, IntersectionArray)

"""
    array(ia::IntersectionArray)

Return the array of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of sets

### Output

The array of an intersection of a finite number of sets.
"""
function array(ia::IntersectionArray)
    return ia.array
end

include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("ispolyhedraltype.jl")
include("in.jl")
include("support_vector.jl")
include("translate.jl")
