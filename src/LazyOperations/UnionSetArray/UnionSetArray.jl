"""
   UnionSetArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the set union of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The union of convex sets is typically not convex.
"""
struct UnionSetArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    UnionSet!(X, Y)

Convenience function to compute the lazy union and modify `UnionSetArray`s
in-place.
"""
function UnionSet! end

# constructor for an empty union with optional size hint and numeric type
function UnionSetArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return UnionSetArray(arr)
end

# EmptySet is the neutral element for UnionSetArray
@neutral(UnionSetArray, EmptySet)

# Universe is the absorbing element for UnionSetArray
@absorbing(UnionSetArray, Universe)

# add functions connecting UnionSet and UnionSetArray
@declare_array_version(UnionSet, UnionSetArray)

"""
   array(cup::UnionSetArray)

Return the array of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The array of the union.
"""
function array(cup::UnionSetArray)
    return cup.array
end

include("an_element.jl")
include("concretize.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
