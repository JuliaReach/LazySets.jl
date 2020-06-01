import Base.isempty

export ConvexHullArray, CHArray,
       array

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
isconvextype(::Type{<:ConvexHullArray}) = true

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

"""
    CHArray

Alias for `ConvexHullArray`.
"""
const CHArray = ConvexHullArray

"""
    array(cha::ConvexHullArray{N, S}) where {N<:Real, S<:LazySet{N}}

Return the array of a convex hull of a finite number of convex sets.

### Input

- `cha` -- convex hull array

### Output

The array of a convex hull of a finite number of convex sets.
"""
function array(cha::ConvexHullArray{N, S}) where {N<:Real, S<:LazySet{N}}
    return cha.array
end

"""
    dim(cha::ConvexHullArray)

Return the dimension of the convex hull of a finite number of convex sets.

### Input

- `cha` -- convex hull array

### Output

The ambient dimension of the convex hull of a finite number of convex sets.
"""
function dim(cha::ConvexHullArray)
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
    isbounded(cha::ConvexHullArray)

Determine whether a convex hull of a finite number of convex sets is
bounded.

### Input

- `cha` -- convex hull of a finite number of convex sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cha::ConvexHullArray)
    return all(isbounded, cha.array)
end

"""
    isempty(cha::ConvexHullArray)

Return if a convex hull array is empty or not.

### Input

- `cha` -- convex hull array

### Output

`true` iff all wrapped sets are empty.
"""
function isempty(cha::ConvexHullArray)
    return all(isempty, array(cha))
end

"""
    vertices_list(cha::ConvexHullArray; apply_convex_hull::Bool=true,
                  backend=nothing)

Return the list of vertices of the convex hull of a finite number of convex sets.

### Input

- `cha`               -- convex hull of a finite number of convex sets
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

The list of vertices.
"""
function vertices_list(cha::ConvexHullArray;
                       apply_convex_hull::Bool=true,
                       backend=nothing)
    vlist = vcat([vertices_list(Xi) for Xi in array(cha)]...)
    if apply_convex_hull
        convex_hull!(vlist, backend=backend)
    end
    return vlist
end

# list of constraints of the convex hull array of singletons
function constraints_list(X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    n = dim(X)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return constraints_list(V)
end

# membership in convex hull array of singletons
function ∈(x::AbstractVector{N}, X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    n = length(x)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return x ∈ V
end
