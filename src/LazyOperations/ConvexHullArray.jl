import Base.isempty

export ConvexHullArray, CHArray,
       array

# ================================
#  Convex hull of an array of sets
# ================================

"""
    ConvexHullArray{N, S<:ConvexSet{N}} <: ConvexSet{N}

Type that represents the symbolic convex hull of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the neutral element for `ConvexHullArray`.

Constructors:

- `ConvexHullArray(array::Vector{<:ConvexSet})` -- default constructor

- `ConvexHullArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty hull with optional size hint and numeric type

### Examples

Convex hull of 100 two-dimensional balls whose centers follow a sinusoidal:

```jldoctest
julia> b = [Ball2([2*pi*i/100, sin(2*pi*i/100)], 0.05) for i in 1:100];

julia> c = ConvexHullArray(b);
```
"""
struct ConvexHullArray{N, S<:ConvexSet{N}} <: ConvexSet{N}
    array::Vector{S}
end

isoperationtype(::Type{<:ConvexHullArray}) = true
isconvextype(::Type{<:ConvexHullArray}) = true

# constructor for an empty hull with optional size hint and numeric type
function ConvexHullArray(n::Int=0, N::Type=Float64)
    a = Vector{ConvexSet{N}}()
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
    array(cha::ConvexHullArray)

Return the array of a convex hull of a finite number of sets.

### Input

- `cha` -- convex hull array

### Output

The array of a convex hull of a finite number of sets.
"""
function array(cha::ConvexHullArray)
    return cha.array
end

"""
    dim(cha::ConvexHullArray)

Return the dimension of the convex hull of a finite number of sets.

### Input

- `cha` -- convex hull array

### Output

The ambient dimension of the convex hull of a finite number of sets.
"""
function dim(cha::ConvexHullArray)
    @assert !isempty(cha.array) "an empty convex hull is not allowed"
    return dim(cha.array[1])
end

"""
    σ(d::AbstractVector, cha::ConvexHullArray)

Return the support vector of a convex hull array in a given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull array
"""
function σ(d::AbstractVector, cha::ConvexHullArray)
    @assert !isempty(cha.array) "an empty convex hull is not allowed"
    svec = d
    N = eltype(d)
    rmax = N(-Inf)
    for chi in cha.array
        si = σ(d, chi)
        ri = dot(d, si)
        if ri > rmax
            rmax = ri
            svec = si
        end
    end
    return svec
end

"""
    ρ(d::AbstractVector, cha::ConvexHullArray)

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
function ρ(d::AbstractVector, cha::ConvexHullArray)
    return maximum(ρ(d, Xi) for Xi in array(cha))
end

"""
    isbounded(cha::ConvexHullArray)

Determine whether a convex hull of a finite number of sets is bounded.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cha::ConvexHullArray)
    return all(isbounded, cha.array)
end

function isboundedtype(::Type{<:ConvexHullArray{N, S}}) where {N, S}
    return isboundedtype(S)
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

Return the list of vertices of the convex hull of a finite number of sets.

### Input

- `cha`               -- convex hull of a finite number of sets
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
function ∈(x::AbstractVector, X::ConvexHullArray)
    n = length(x)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return x ∈ V
end

function concretize(cha::ConvexHullArray)
    a = array(cha)
    @assert !isempty(a) "an empty convex hull is not allowed"
    X = cha
    @inbounds for (i, Y) in enumerate(a)
        if i == 1
            X = concretize(Y)
        else
            X = convex_hull(X, concretize(Y))
        end
    end
    return X
end
