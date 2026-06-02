"""
   MinkowskiSumArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSumArray`.

The Minkowski sum preserves convexity: if the set arguments are convex, then
their Minkowski sum is convex as well.

The convenience aliases `⊕` and `+` are also available. `⊕` can be typed by
`\\oplus<tab>`.
"""
struct MinkowskiSumArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    MinkowskiSum!(X, Y)

Convenience function to compute the lazy Minkowski sum and modify
`MinkowskiSumArray`s in-place.
"""
function MinkowskiSum! end

+(X::LazySet, Xs::LazySet...) = MinkowskiSumArray(vcat(X, Xs...))
+(X::LazySet) = X
+(Xs::Vector{<:LazySet}) = MinkowskiSumArray(Xs)

⊕(X::LazySet, Xs::LazySet...) = +(X, Xs...)
⊕(Xs::Vector{<:LazySet}) = +(Xs)

# constructor for an empty sum with optional size hint and numeric type
function MinkowskiSumArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

# ZeroSet is the neutral element for MinkowskiSumArray
@neutral(MinkowskiSumArray, ZeroSet)

# EmptySet and Universe are the absorbing elements for MinkowskiSumArray
@absorbing(MinkowskiSumArray, EmptySet)
# @absorbing(MinkowskiSumArray, Universe)  # TODO problematic

# add functions connecting MinkowskiSum and MinkowskiSumArray
@declare_array_version(MinkowskiSum, MinkowskiSumArray)

"""
   array(msa::MinkowskiSumArray)

Return the array of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

The array of a Minkowski sum of a finite number of sets.
"""
function array(msa::MinkowskiSumArray)
    return msa.array
end

include("center.jl")
include("concretize.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
