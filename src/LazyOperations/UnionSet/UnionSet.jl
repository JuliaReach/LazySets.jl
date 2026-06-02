import Base: ∪

"""
    UnionSet{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the set union of two sets.

### Fields

- `X` -- set
- `Y` -- set

### Notes

The union of convex sets is typically not convex.

The convenience alias `∪` can be typed by `\\cup<tab>`.
"""
struct UnionSet{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function UnionSet(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a union must have the same dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

# EmptySet is the neutral element for UnionSet
@neutral(UnionSet, EmptySet)

# Universe is the absorbing element for UnionSet
@absorbing(UnionSet, Universe)

# interface for binary set operations
first(U::UnionSet) = U.X
second(U::UnionSet) = U.Y
@declare_binary_operation(UnionSet)

∪(X::LazySet, Y::LazySet) = UnionSet(X, Y)

"""
    swap(cup::UnionSet)

Return a new `UnionSet` object with the arguments swapped.

### Input

- `cup` -- union of two sets

### Output

A new `UnionSet` object with the arguments swapped.
"""
function swap(cup::UnionSet)
    return UnionSet(cup.Y, cup.X)
end

include("an_element.jl")
include("concretize.jl")
include("convex_hull.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("vertices_list.jl")
include("volume.jl")
include("in.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
