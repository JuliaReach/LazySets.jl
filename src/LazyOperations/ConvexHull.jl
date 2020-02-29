import Base.isempty

export ConvexHull, CH,
       convex_hull,
       convex_hull!,
       ConvexHull!,
       swap

"""
    ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the convex hull of the union of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set

### Notes

The `EmptySet` is the neutral element for `ConvexHull`.

### Examples

Convex hull of two 100-dimensional Euclidean balls:

```jldoctest
julia> b1, b2 = Ball2(zeros(100), 0.1), Ball2(4*ones(100), 0.2);

julia> c = ConvexHull(b1, b2);

julia> typeof(c)
ConvexHull{Float64,Ball2{Float64,Array{Float64,1}},Ball2{Float64,Array{Float64,1}}}
```
"""
struct ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ConvexHull(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N},
                                             S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in a convex hull must have the same " *
            "dimension"
        return new{N, S1, S2}(X, Y)
    end
end

isoperationtype(::Type{<:ConvexHull}) = true
isconvextype(::Type{<:ConvexHull}) = true

# EmptySet is the neutral element for ConvexHull
@neutral(ConvexHull, EmptySet)

# Universe is the absorbing element for ConvexHull
@absorbing(ConvexHull, Universe)

"""
    CH

Alias for `ConvexHull`.
"""
const CH = ConvexHull

"""
    swap(ch::ConvexHull)

Return a new `ConvexHull` object with the arguments swapped.

### Input

- `ch` -- convex hull of two convex sets

### Output

A new `ConvexHull` object with the arguments swapped.
"""
function swap(ch::ConvexHull)
    return ConvexHull(ch.Y, ch.X)
end

"""
    dim(ch::ConvexHull)

Return the dimension of a convex hull of two convex sets.

### Input

- `ch` -- convex hull of two convex sets

### Output

The ambient dimension of the convex hull of two convex sets.
"""
function dim(ch::ConvexHull)
    return dim(ch.X)
end

"""
    σ(d::AbstractVector{N}, ch::ConvexHull{N}) where {N<:Real}

Return the support vector of a convex hull of two convex sets in a given
direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two convex sets

### Output

The support vector of the convex hull in the given direction.
"""
function σ(d::AbstractVector{N}, ch::ConvexHull{N}) where {N<:Real}
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
    ρ1 = dot(d, σ1)
    ρ2 = dot(d, σ2)
    return ρ1 >= ρ2 ? σ1 : σ2
end

"""
    ρ(d::AbstractVector{N}, ch::ConvexHull{N}) where {N<:Real}

Return the support function of a convex hull of two convex sets in a given
direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two convex sets

### Output

The support function of the convex hull in the given direction.
"""
function ρ(d::AbstractVector{N}, ch::ConvexHull{N}) where {N<:Real}
    return max(ρ(d, ch.X), ρ(d, ch.Y))
end

"""
    isbounded(ch::ConvexHull)

Determine whether a convex hull of two convex sets is bounded.

### Input

- `ch` -- convex hull of two convex sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ch::ConvexHull)
    return isbounded(ch.X) && isbounded(ch.Y)
end

"""
    isempty(ch::ConvexHull)

Return if a convex hull of two convex sets is empty or not.

### Input

- `ch` -- convex hull

### Output

`true` iff both wrapped sets are empty.
"""
function isempty(ch::ConvexHull)
    return isempty(ch.X) && isempty(ch.Y)
end
