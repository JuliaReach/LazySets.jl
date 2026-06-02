import Base: +, getindex

"""
    MinkowskiSum{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of two sets, i.e., the set

```math
X ⊕ Y = \\{x + y : x ∈ X, y ∈ Y\\}.
```

### Fields

- `X` -- set
- `Y` -- set

### Notes

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSum`.

The Minkowski sum preserves convexity: if the set arguments are convex, then
their Minkowski sum is convex as well.

The convenience aliases `⊕` and `+` are also available. `⊕` can be typed by
`\\oplus<tab>`.
"""
struct MinkowskiSum{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function MinkowskiSum(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a Minkowski sum must have the same dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

+(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

⊕(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

concrete_function(::Type{<:MinkowskiSum}) = minkowski_sum

# ZeroSet is the neutral element for MinkowskiSum
@neutral(MinkowskiSum, ZeroSet)

# EmptySet and Universe are the absorbing elements for MinkowskiSum
@absorbing(MinkowskiSum, EmptySet)
# @absorbing(MinkowskiSum, Universe)  # TODO problematic

# interface for binary set operations
first(ms::MinkowskiSum) = ms.X
second(ms::MinkowskiSum) = ms.Y
@declare_binary_operation(MinkowskiSum)

"""
    swap(ms::MinkowskiSum)

Return a new `MinkowskiSum` object with the arguments swapped.

### Input

- `ms` -- Minkowski sum of two sets

### Output

A new `MinkowskiSum` object with the arguments swapped.
"""
function swap(ms::MinkowskiSum)
    return MinkowskiSum(ms.Y, ms.X)
end

include("center.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("ispolyhedraltype.jl")
include("vertices_list.jl")
include("in.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
