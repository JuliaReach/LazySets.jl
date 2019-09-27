import Base.isempty

export ConvexHull, CH,
       convex_hull,
       convex_hull!,
       ConvexHullArray, CHArray,
       ConvexHull!,
       array,
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
ConvexHull{Float64,Ball2{Float64},Ball2{Float64}}
```
"""
struct ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ConvexHull{N, S1, S2}(X::S1, Y::S2) where
            {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in a convex hull must have the same " *
            "dimension"
        return new{N, S1, S2}(X, Y)
    end
end

isoperationtype(::Type{<:ConvexHull}) = true

# convenience constructor without type parameter
ConvexHull(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} =
    ConvexHull{N, S1, S2}(X, Y)

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
    dim(ch::ConvexHull)::Int

Return the dimension of a convex hull of two convex sets.

### Input

- `ch` -- convex hull of two convex sets

### Output

The ambient dimension of the convex hull of two convex sets.
"""
function dim(ch::ConvexHull)::Int
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
    isbounded(ch::ConvexHull)::Bool

Determine whether a convex hull of two convex sets is bounded.

### Input

- `ch` -- convex hull of two convex sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ch::ConvexHull)::Bool
    return isbounded(ch.X) && isbounded(ch.Y)
end

"""
    isempty(ch::ConvexHull)::Bool

Return if a convex hull of two convex sets is empty or not.

### Input

- `ch` -- convex hull

### Output

`true` iff both wrapped sets are empty.
"""
function isempty(ch::ConvexHull)::Bool
    return isempty(ch.X) && isempty(ch.Y)
end

# ================================
#  Convex hull of an array of sets
# ================================
"""
    ConvexHullArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the symbolic convex hull of a finite number of convex sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the neutral element for `ConvexHullArray`.

Constructors:

- `ConvexHullArray(array::Vector{<:LazySet})` -- default constructor

- `ConvexHullArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty hull with optional size hint and numeric type

### Examples

Convex hull of 100 two-dimensional balls whose centers follows a sinusoidal:

```jldoctest
julia> b = [Ball2([2*pi*i/100, sin(2*pi*i/100)], 0.05) for i in 1:100];

julia> c = ConvexHullArray(b);
```
"""
struct ConvexHullArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

isoperationtype(::Type{<:ConvexHullArray}) = true

# constructor for an empty hull with optional size hint and numeric type
function ConvexHullArray(n::Int=0, N::Type=Float64)::ConvexHullArray
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

"""
    CHArray

Alias for `ConvexHullArray`.
"""
const CHArray = ConvexHullArray

"""
    array(cha::ConvexHullArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a convex hull of a finite number of convex sets.

### Input

- `cha` -- convex hull array

### Output

The array of a convex hull of a finite number of convex sets.
"""
function array(cha::ConvexHullArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cha.array
end

"""
    dim(cha::ConvexHullArray)::Int

Return the dimension of the convex hull of a finite number of convex sets.

### Input

- `cha` -- convex hull array

### Output

The ambient dimension of the convex hull of a finite number of convex sets.
"""
function dim(cha::ConvexHullArray)::Int
    @assert !isempty(cha.array)
    return dim(cha.array[1])
end

"""
    σ(d::AbstractVector{N}, cha::ConvexHullArray{N}) where {N<:Real}

Return the support vector of a convex hull array in a given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull array
"""
function σ(d::AbstractVector{N}, cha::ConvexHullArray{N}) where {N<:Real}
    s = σ(d, cha.array[1])
    ri = dot(d, s)
    rmax = ri
    for (i, chi) in enumerate(cha.array[2:end])
        si = σ(d, chi)
        ri = dot(d, si)
        if ri > rmax
            rmax = ri
            s = si
        end
    end
    return s
end

"""
    ρ(d::AbstractVector{N}, cha::ConvexHullArray{N}) where {N<:Real}

Return the support function of a convex hull array in a given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull array

### Output

The support function of the convex hull array in the given direction.

### Algorithm

This algorihm calculates the maximum over all ``ρ(d, X_i)`` where the
``X_1, …, X_k`` are the sets in the array `cha`.
"""
function ρ(d::AbstractVector{N}, cha::ConvexHullArray{N}) where {N<:Real}
    return maximum([ρ(d, Xi) for Xi in array(cha)])
end

"""
    isbounded(cha::ConvexHullArray)::Bool

Determine whether a convex hull of a finite number of convex sets is
bounded.

### Input

- `cha` -- convex hull of a finite number of convex sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cha::ConvexHullArray)::Bool
    return all(x -> isbounded(x), cha.array)
end

"""
    isempty(cha::ConvexHullArray)::Bool

Return if a convex hull array is empty or not.

### Input

- `cha` -- convex hull array

### Output

`true` iff all wrapped sets are empty.
"""
function isempty(cha::ConvexHullArray)::Bool
    return all(X -> isempty(X), array(cha))
end
