import Base.<=

export AbstractPolygon,
       tohrep,
       tovrep

"""
    AbstractPolygon{N<:Real} <: AbstractPolytope{N}

Abstract type for polygons (i.e., 2D polytopes).

### Notes

Every concrete `AbstractPolygon` must define the following functions:
- `tovrep(::AbstractPolygon{N})::VPolygon{N}`         -- transform into
    V-representation
- `tohrep(::AbstractPolygon{N})::AbstractHPolygon{N}` -- transform into
    H-representation

```jldoctest
julia> subtypes(AbstractPolygon)
2-element Array{Any,1}:
 AbstractHPolygon
 VPolygon
```
"""
abstract type AbstractPolygon{N<:Real} <: AbstractPolytope{N} end


# --- LazySet interface functions ---


"""
    dim(P::AbstractPolygon)::Int

Return the ambient dimension of a polygon.

### Input

- `P` -- polygon

### Output

The ambient dimension of the polygon, which is 2.
"""
@inline function dim(P::AbstractPolygon)::Int
    return 2
end

"""
    jump2pi(x::N)::N where {N<:AbstractFloat}

Return ``x + 2π`` if ``x`` is negative, otherwise return ``x``.

### Input

- `x` -- real scalar

### Output

``x + 2π`` if ``x`` is negative, ``x`` otherwise.

### Examples

```jldoctest
julia> import LazySets.jump2pi

julia> jump2pi(0.0)
0.0

julia> jump2pi(-0.5)
5.783185307179586

julia> jump2pi(0.5)
0.5
```
"""
@inline function jump2pi(x::N)::N where {N<:AbstractFloat}
    x < zero(N) ? 2 * pi + x : x
end

"""
    quadrant(w::AbstractVector{N})::Int where {N<:Real}

Compute the quadrant where the direction `w` belongs.

### Input

- `w` --  direction

### Output

An integer from 0 to 3, with the following convention:

```
     ^
   1 | 0
  ---+-->
   2 | 3
```

### Algorithm

The idea is to encode the following logic function:
``11 ↦ 0, 01 ↦ 1, 00 ↦ 2, 10 ↦ 3``, according to the convention of above.

This function is inspired from AGPX's answer in:
[Sort points in clockwise order?](https://stackoverflow.com/a/46635372)
"""
@inline function quadrant(w::AbstractVector{N})::Int where {N<:Real}
    dwx = w[1] >= zero(N) ? 1 : 0
    dwy = w[2] >= zero(N) ? 1 : 0
    return (1 - dwx) + (1 - dwy) + ((dwx & (1 - dwy)) << 1)
end

"""
    <=(u::AbstractVector{N}, v::AbstractVector{N})::Bool where {N<:Real}

Compares two 2D vectors by their direction.

### Input

- `u` --  first 2D direction
- `v` --  second 2D direction

### Output

True iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).

### Algorithm

The implementation checks the quadrant of each direction, and compares directions
using the right-hand rule. In particular, it doesn't use the arctangent.
"""
function <=(u::AbstractVector{N}, v::AbstractVector{N})::Bool where {N<:Real}
    qu, qv = quadrant(u), quadrant(v)
    if qu == qv
        # same quadrant, check right-turn with center 0
        ans = u[1] * v[2] - v[1] * u[2] >= zero(N)
    else
        # different quadrant
        ans = qu < qv
    end
    return ans
end

"""
    <=(u::AbstractVector{N}, v::AbstractVector{N})::Bool where {N<:AbstractFloat}

Compares two 2D vectors by their direction.

### Input

- `u` --  first 2D direction
- `v` --  second 2D direction

### Output

True iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).

### Algorithm

The implementation uses the arctangent function with sign, `atan`, which for two
arguments implements the
[`atan2` function](https://en.wikipedia.org/wiki/Atan2).
"""
function <=(u::AbstractVector{N},
            v::AbstractVector{N})::Bool where {N<:AbstractFloat}
    return jump2pi(atan(u[2], u[1])) <= jump2pi(atan(v[2], v[1]))
end

"""
    linear_map(M::AbstractMatrix, P::AbstractPolygon{N};
               output_type::Type{<:LazySet}=typeof(P)) where {N}

Concrete linear map of an abstract polygon.

### Input

- `M`           -- matrix
- `P`           -- abstract polygon
- `output_type` -- (optional, default: type of `P`) type of the result

### Output

A set of type `output_type`.

### Algorithm

The linear map ``M`` is applied to each vertex of the given set ``P``,
obtaining a polygon in V-representation. Since polygons are closed under linear
map, by default ``MP`` is converted to the concrete type of ``P``. If an
`output_type` is given, the corresponding `convert` method is invoked.
"""
function linear_map(M::AbstractMatrix, P::AbstractPolygon{N};
                    output_type::Type{<:LazySet}=typeof(P)) where {N}
    @assert dim(P) == size(M, 2)
    MP = broadcast(v -> M * v, vertices_list(P)) |> VPolygon{N}
    return convert(output_type, MP)
end
