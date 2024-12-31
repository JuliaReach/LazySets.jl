export ConvexHullArray, CHArray,
       array

"""
    ConvexHullArray{N, S<:LazySet{N}} <: ConvexSet{N}

Type that represents the symbolic convex hull of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the neutral element for `ConvexHullArray`.

A `ConvexHullArray` is always convex.

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

isoperationtype(::Type{<:ConvexHullArray}) = true

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

"""
    dim(cha::ConvexHullArray)

Return the dimension of the convex hull of a finite number of sets.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

The ambient dimension of the convex hull of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(cha::ConvexHullArray)
    return length(cha.array) == 0 ? 0 : dim(cha.array[1])
end

"""
    σ(d::AbstractVector, cha::ConvexHullArray)

Return a support vector of a convex hull of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `cha` -- convex hull of a finite number of sets

### Output

A support vector in the given direction.
"""
function σ(d::AbstractVector, cha::ConvexHullArray)
    @assert !isempty(cha.array) "an empty convex hull is not allowed"
    return _σ_union(d, array(cha))
end

"""
    ρ(d::AbstractVector, cha::ConvexHullArray)

Evaluate the support function of a convex hull of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull of a finite number of sets

### Output

The evaluation of the support function of the convex hull of a finite number of
sets in the given direction.

### Algorithm

This algorithm calculates the maximum over all ``ρ(d, X_i)``, where the
``X_1, …, X_k`` are the sets in the array of `cha`.
"""
function ρ(d::AbstractVector, cha::ConvexHullArray)
    return maximum(ρ(d, Xi) for Xi in cha)
end

"""
    isbounded(cha::ConvexHullArray)

Check whether a convex hull of a finite number of sets is bounded.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cha::ConvexHullArray)
    return all(isbounded, cha.array)
end

function isboundedtype(::Type{<:ConvexHullArray{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
    isempty(cha::ConvexHullArray)

Check whether a convex hull of a finite number of sets is empty.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

`true` iff all wrapped sets are empty.
"""
function isempty(cha::ConvexHullArray)
    return all(isempty, array(cha))
end

"""
    vertices_list(cha::ConvexHullArray; [apply_convex_hull]::Bool=true,
                  [backend]=nothing, [prune]::Bool=apply_convex_hull)

Return a list of vertices of the convex hull of a finite number of sets.

### Input

- `cha`               -- convex hull of a finite number of sets
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)
- `prune`             -- (optional, default: `apply_convex_hull`) alias for
                         `apply_convex_hull`

### Output

A list of vertices.
"""
function vertices_list(cha::ConvexHullArray;
                       apply_convex_hull::Bool=true,
                       backend=nothing,
                       prune::Bool=apply_convex_hull)
    vlist = vcat([vertices_list(Xi) for Xi in cha]...)
    if apply_convex_hull || prune
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end

# list of constraints of a convex-hull array of singletons
function constraints_list(X::ConvexHullArray{N,Singleton{N,VT}}) where {N,VT}
    n = dim(X)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return constraints_list(V)
end

# membership in a convex-hull array of singletons
function ∈(x::AbstractVector, X::ConvexHullArray)
    n = length(x)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return x ∈ V
end

function concretize(cha::ConvexHullArray)
    a = array(cha)
    @assert !isempty(a) "an empty convex hull is not allowed"
    if all(ispolyhedral, a)
        return _convex_hull_polytopes(cha)
    end

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

function translate(cha::ConvexHullArray, x::AbstractVector)
    return ConvexHullArray([translate(X, x) for X in array(cha)])
end
