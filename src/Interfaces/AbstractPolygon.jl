export AbstractPolygon,
       tohrep,
       tovrep

"""
    AbstractPolygon{N} <: AbstractPolytope{N}

Abstract type for convex polygons (i.e., two-dimensional polytopes).

### Notes

Every concrete `AbstractPolygon` must define the following functions:

- `tohrep(::AbstractPolygon)` -- transform into constraint representation
- `tovrep(::AbstractPolygon)` -- transform into vertex representation

The subtypes of `AbstractPolygon` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolygon)
2-element Vector{Any}:
 AbstractHPolygon
 VPolygon
```
"""
abstract type AbstractPolygon{N} <: AbstractPolytope{N} end

"""
    tohrep(P::AbstractPolygon)

Convert a convex polygon to constraint representation.

### Input

- `P` -- convex polygon

### Output

A polygon in constraint representation.
"""
function tohrep(::AbstractPolygon) end

"""
    tovrep(P::AbstractPolygon)

Convert a convex polygon to vertex representation.

### Input

- `P` -- convex polygon

### Output

A polygon in vertex representation.
"""
function tovrep(::AbstractPolygon) end

isconvextype(::Type{<:AbstractPolygon}) = true

"""
    dim(P::AbstractPolygon)

Return the ambient dimension of a convex polygon.

### Input

- `P` -- convex polygon

### Output

The ambient dimension of the polygon, which is 2.
"""
@inline function dim(::AbstractPolygon)
    return 2
end

"""
    jump2pi(x::N) where {N<:AbstractFloat}

Return ``x + 2π`` if ``x`` is negative, otherwise return ``x``.

### Input

- `x` -- real scalar

### Output

``x + 2π`` if ``x`` is negative, ``x`` otherwise.

### Examples

```jldoctest
julia> using LazySets: jump2pi

julia> jump2pi(0.0)
0.0

julia> jump2pi(-0.5)
5.783185307179586

julia> jump2pi(0.5)
0.5
```
"""
@inline function jump2pi(x::N) where {N<:AbstractFloat}
    return x < zero(N) ? 2 * pi + x : x
end

"""
    quadrant(w::AbstractVector{N}) where {N}

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
``11 ↦ 0, 01 ↦ 1, 00 ↦ 2, 10 ↦ 3``, according to the convention above.

This function is inspired from AGPX's answer in:
[Sort points in clockwise order?](https://stackoverflow.com/a/46635372)
"""
@inline function quadrant(w::AbstractVector{N}) where {N}
    dwx = w[1] >= zero(N) ? 1 : 0
    dwy = w[2] >= zero(N) ? 1 : 0
    return (1 - dwx) + (1 - dwy) + ((dwx & (1 - dwy)) << 1)
end

"""
    ⪯(u::AbstractVector, v::AbstractVector)

Compare two 2D vectors by their direction.

### Input

- `u` -- first 2D direction
- `v` -- second 2D direction

### Output

`true` iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``.

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).

### Algorithm

The implementation checks the quadrant of each direction, and compares
directions using the right-hand rule.
In particular, this method does not use the arctangent.
"""
function ⪯(u::AbstractVector, v::AbstractVector)
    @assert length(u) == length(v) == 2 "comparison of vectors `u` and `v` " *
                                        "by their direction requires that they are of length 2, " *
                                        "but their lengths are $(length(u)) and $(length(v)), respectively"

    return _leq_quadrant(u, v)
end

function _leq_quadrant(u::AbstractVector, v::AbstractVector)
    qu, qv = quadrant(u), quadrant(v)
    if qu == qv
        # same quadrant, check right-turn with center 0
        return is_right_turn(u, v)
    end
    # different quadrant
    return qu < qv
end

"""
    _leq_trig(u::AbstractVector{N}, v::AbstractVector{N}) where {N<:AbstractFloat}

Compare two 2D vectors by their direction.

### Input

- `u` --  first 2D direction
- `v` --  second 2D direction

### Output

`true` iff ``\\arg(u) [2π] ≤ \\arg(v) [2π]``.

### Notes

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).

### Algorithm

The implementation uses the arctangent function with sign, `atan`, which for two
arguments implements the [`atan2` function](https://en.wikipedia.org/wiki/Atan2).
"""
function _leq_trig(u::AbstractVector{N}, v::AbstractVector{N}) where {N<:AbstractFloat}
    return jump2pi(atan(u[2], u[1])) <= jump2pi(atan(v[2], v[1]))
end

"""
    volume(P::AbstractPolygon)

Compute the volume of a convex polygon.

### Input

- `P` -- convex polygon

### Output

A number representing the volume of `P`.

### Notes

In 2D the volume is equivalent to the area.
"""
function volume(P::AbstractPolygon)
    return area(P)
end

# Note: this method assumes that the vertices are sorted in CCW order
function _intersection_vrep_2d(spoly::AbstractVector{VT},
                               cpoly::AbstractVector{VT}) where {VT<:AbstractVector}
    outarr = spoly
    q = cpoly[end]
    for p in cpoly
        inarr = outarr
        outarr = Vector{VT}()
        isempty(inarr) && break
        s = inarr[end]
        for vtx in inarr
            if _isinside(vtx, q, p)
                if !_isinside(s, q, p)
                    push!(outarr, _intersection_line_segments(q, p, s, vtx))
                end
                push!(outarr, vtx)
            elseif _isinside(s, q, p)
                push!(outarr, _intersection_line_segments(q, p, s, vtx))
            end
            s = vtx
        end
        q = p
    end
    return outarr
end

function _isinside(p, a, b)
    @inbounds begin
        α = (b[1] - a[1]) * (p[2] - a[2])
        β = (b[2] - a[2]) * (p[1] - a[1])
    end
    return !_leq(α, β)
end

function _intersection_line_segments(a, b, s, f)
    @inbounds begin
        dc = [a[1] - b[1], a[2] - b[2]]
        dp = [s[1] - f[1], s[2] - f[2]]
        n1 = a[1] * b[2] - a[2] * b[1]
        n2 = s[1] * f[2] - s[2] * f[1]
        n3 = 1.0 / (dc[1] * dp[2] - dc[2] * dp[1])

        α = (n1 * dp[1] - n2 * dc[1]) * n3
        β = (n1 * dp[2] - n2 * dc[2]) * n3
    end
    return [α, β]
end
