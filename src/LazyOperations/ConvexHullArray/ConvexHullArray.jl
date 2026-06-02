"""
    ConvexHullArray{N, S<:LazySet{N}} <: ConvexSet{N}

Type that represents the symbolic convex hull of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the neutral element for `ConvexHullArray`.

A `ConvexHullArray` is always convex.

The convenience alias `CHArray` is also available.

### Examples

Convex hull of 100 two-dimensional balls whose centers follow a sinusoidal:

```jldoctest
julia> b = [Ball2([2*pi*i/100, sin(2*pi*i/100)], 0.05) for i in 1:100];

julia> c = ConvexHullArray(b);
```
"""
struct ConvexHullArray{N,S<:LazySet{N}} <: ConvexSet{N}
    array::Vector{S}
end

"""
    ConvexHull!(X, Y)

Convenience function to compute the lazy convex hull and modify
`ConvexHullArray`s in-place.
"""
function ConvexHull! end

# constructor for an empty hull with optional size hint and numeric type
function ConvexHullArray(n::Int=0, N::Type=Float64)
    a = Vector{LazySet{N}}()
    sizehint!(a, n)
    return ConvexHullArray(a)
end

# EmptySet is the neutral element for ConvexHullArray
@neutral(ConvexHullArray, EmptySet)

# Universe is the absorbing element for ConvexHullArray
@absorbing(ConvexHullArray, Universe)

# add functions connecting ConvexHull and ConvexHullArray
@declare_array_version(ConvexHull, ConvexHullArray)

const CHArray = ConvexHullArray

"""
    array(cha::ConvexHullArray)

Return the array of a convex hull of a finite number of sets.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

The array of a convex hull of a finite number of sets.
"""
function array(cha::ConvexHullArray)
    return cha.array
end

include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("ispolyhedraltype.jl")
include("vertices_list.jl")
include("in.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
