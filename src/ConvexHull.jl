include("convex_hull_algorithms.jl")

export ConvexHull, CH,
       convex_hull,
       convex_hull!,
       ConvexHullArray, CHArray,
       array

"""
    ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the convex hull of the union of two sets.

### Fields

- `X` -- first set
- `Y` -- second set

### Notes

The `EmptySet` is the neutral element for `ConvexHull`.

### Examples

Convex hull of two 100-dimensional Euclidean balls:

```jldoctest
julia> b1, b2 = Ball2(zeros(100), 0.1), Ball2(4*ones(100), 0.2);

julia> c = ConvexHull(b1, b2);

julia> typeof(c)
LazySets.ConvexHull{Float64,LazySets.Ball2{Float64},LazySets.Ball2{Float64}}
```
"""
struct ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    ConvexHull{N, S1, S2}(X::S1, Y::S2) where
        {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
            dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
ConvexHull(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    ConvexHull{N, S1, S2}(X, Y)

# EmptySet is the neutral element for ConvexHull
@neutral(ConvexHull, EmptySet)

"""
    CH

Alias for `ConvexHull`.
"""
const CH = ConvexHull

"""
    dim(ch::ConvexHull)::Int

Return the dimension of a convex hull of two sets.

### Input

- `ch` -- convex hull of two sets

### Output

The ambient dimension of the convex hull of two sets.
"""
function dim(ch::ConvexHull)::Int
    return dim(ch.X)
end

"""
    σ(d::AbstractVector{<:Real}, ch::ConvexHull)::AbstractVector{<:Real}

Return the support vector of a convex hull of two sets in a given direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two sets
"""
function σ(d::AbstractVector{<:Real}, ch::ConvexHull)::AbstractVector{<:Real}
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
    ρ1 = dot(d, σ1)
    ρ2 = dot(d, σ2)
    return ρ1 >= ρ2 ? σ1 : σ2
end

# ================================
#  Convex hull of an array of sets
# ================================
"""
    ConvexHullArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the symbolic convex hull of a finite number of sets.

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

# type-less convenience constructor
ConvexHullArray(a::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
    ConvexHullArray{N, S}(a)

# constructor for an empty hull with optional size hint and numeric type
function ConvexHullArray(n::Int=0, N::Type=Float64)::ConvexHullArray
    a = Vector{LazySet{N}}(0)
    sizehint!(a, n)
    return ConvexHullArray(a)
end

# EmptySet is the neutral element for ConvexHullArray
@neutral(ConvexHullArray, EmptySet)

# add functions connecting ConvexHull and ConvexHullArray
@declare_array_version(ConvexHull, ConvexHullArray)

"""
    CHArray

Alias for `ConvexHullArray`.
"""
const CHArray = ConvexHullArray

"""
    array(cha::ConvexHullArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a convex hull of a finite number of sets.

### Input

- `cha` -- convex hull array

### Output

The array of a convex hull of a finite number of sets.
"""
function array(cha::ConvexHullArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cha.array
end

"""
    dim(cha::ConvexHullArray)::Int

Return the dimension of the convex hull of a finite number of sets.

### Input

- `cha` -- convex hull array

### Output

The ambient dimension of the convex hull of a finite number of sets.
"""
function dim(cha::ConvexHullArray)::Int
    @assert !isempty(cha.array)
    return dim(cha.array[1])
end

"""
    σ(d::AbstractVector{<:Real}, cha::ConvexHullArray)::AbstractVector{<:Real}

Return the support vector of a convex hull array in a given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull array
"""
function σ(d::AbstractVector{<:Real}, cha::ConvexHullArray)::AbstractVector{<:Real}
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
