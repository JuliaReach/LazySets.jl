"""
    ConvexHull{N, S1<:LazySet{N}, S2<:LazySet{N}} <: ConvexSet{N}

Type that represents the convex hull of the union of two sets, i.e., the set

```math
Z = \\{z ∈ ℝ^n : z = λx + (1-λ)y,\\qquad x ∈ X, y ∈ Y,λ ∈ [0, 1] \\}.
```

### Fields

- `X` -- set
- `Y` -- set

### Notes

The `EmptySet` is the neutral element for `ConvexHull`.

This type is always convex.

The convenience alias `CH` is also available.

### Examples

The convex hull of two 100-dimensional Euclidean balls:

```jldoctest
julia> b1, b2 = Ball2(zeros(100), 0.1), Ball2(4*ones(100), 0.2);

julia> c = ConvexHull(b1, b2);

julia> typeof(c)
ConvexHull{Float64, Ball2{Float64, Vector{Float64}}, Ball2{Float64, Vector{Float64}}}
```
"""
struct ConvexHull{N,S1<:LazySet{N},S2<:LazySet{N}} <: ConvexSet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ConvexHull(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a convex hull must have the same " *
                                 "dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

# constructor with a single argument
ConvexHull(X::LazySet) = ConvexHull(X, X)

concrete_function(::Type{<:ConvexHull}) = convex_hull

# EmptySet is the neutral element for ConvexHull
@neutral(ConvexHull, EmptySet)

# Universe is the absorbing element for ConvexHull
@absorbing(ConvexHull, Universe)

# interface for binary set operations
first(ch::ConvexHull) = ch.X
second(ch::ConvexHull) = ch.Y
@declare_binary_operation(ConvexHull)

const CH = ConvexHull

"""
    swap(ch::ConvexHull)

Return a new `ConvexHull` object with the arguments swapped.

### Input

- `ch` -- convex hull of two sets

### Output

A new `ConvexHull` object with the arguments swapped.
"""
function swap(ch::ConvexHull)
    return ConvexHull(ch.Y, ch.X)
end

include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("ispolyhedraltype.jl")
include("vertices_list.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
