import Base.convert

#= conversion between set types =#

"""
    convert(::Type{HPOLYGON1}, P::HPOLYGON2) where
        {HPOLYGON1<:Union{HPolygon, HPolygonOpt}, HPOLYGON2<:AbstractHPolygon}

Convert between polygon types in H-representation.

### Input

- `type` -- target type
- `P`    -- source polygon

### Output

The polygon represented as the target type.

### Notes

We need the `Union` type for `HPOLYGON1` because the target type must be
concrete.
"""
function convert(::Type{HPOLYGON1},
                 P::HPOLYGON2) where {HPOLYGON1<:Union{HPolygon, HPolygonOpt},
                                      HPOLYGON2<:AbstractHPolygon}
    return HPOLYGON1(P.constraints)
end

"""
    convert(::Type{HPolytope}, P::AbstractHPolygon)

Convert from polygon in H-representation to polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polygon

### Output

The polygon represented as 2D polytope.
"""
function convert(::Type{HPolytope}, P::AbstractHPolygon)
    return HPolytope(P.constraints)
end

"""
    convert(::Type{HPOLYGON}, P::HPolytope) where {HPOLYGON<:AbstractHPolygon}

Convert from 2D polytope in H-representation to polygon in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope (must be 2D)

### Output

The 2D polytope represented as polygon.
"""
function convert(::Type{HPOLYGON},
                 P::HPolytope{N}) where {N, HPOLYGON<:AbstractHPolygon}
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    H = HPOLYGON{N}()
    for ci in constraints_list(P)
        # guarantee that the edges are correctly sorted for storage
        addconstraint!(H, ci)
    end
    return H
end

"""
    convert(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}

Converts a hyperrectangular set to a zonotope.

### Input

- `Zonotope`
- `H` -- hyperrectangular set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    return Zonotope{N}(center(H), Diagonal(radius_hyperrectangle(H)))
end

"""
    convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
           ) where {N, HPOLYGON<:AbstractHPolygon}

Convert from singleton to polygon in H-representation.

### Input

- `type` -- target type
- `S`    -- singleton

### Output

A polygon in constraint representation with the minimal number of constraints
(three).
"""
function convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
                ) where {N, HPOLYGON<:AbstractHPolygon}
    constraints_list = Vector{LinearConstraint{N}}(undef, 3)
    o = one(N)
    z = zero(N)
    v = element(S)
    constraints_list[1] = LinearConstraint([o, o], v[1] + v[2])
    constraints_list[2] = LinearConstraint([-o, z], -v[1])
    constraints_list[3] = LinearConstraint([z, -o], -v[2])
    return HPOLYGON{N}(constraints_list)
end

"""
    convert(::Type{HPOLYGON}, L::LineSegment{N}
          ) where {N, HPOLYGON<:AbstractHPolygon}

Convert from line segment to polygon in H-representation.

### Input

- `type` -- target type
- `L`    -- line segment

### Output

A flat polygon in constraint representation with the minimal number of
constraints (four).
"""
function convert(::Type{HPOLYGON}, L::LineSegment{N}
                ) where {N, HPOLYGON<:AbstractHPolygon}
    H = HPOLYGON{N}()
    c = halfspace_left(L.p, L.q)
    addconstraint!(H, c)
    addconstraint!(H, LinearConstraint(-c.a, -c.b))
    line_dir = L.q - L.p
    c = LinearConstraint(line_dir, dot(L.q, line_dir))
    addconstraint!(H, c)
    line_dir = -line_dir
    addconstraint!(H, LinearConstraint(line_dir, dot(L.p, line_dir)))
    return H
end

import IntervalArithmetic.AbstractInterval

"""
    convert(::Type{Hyperrectangle}, x::Interval{N, IN}) where {N, IN <: AbstractInterval{N}}

Converts a unidimensional interval into a hyperrectangular set.

### Input

- `AbstractHyperrectangle`
- `x` -- interval

### Output

A hyperrectangle.

### Examples

```jldoctest convert_hyperrectangle_interval
julia> convert(Hyperrectangle, Interval(0.0, 1.0))
Hyperrectangle{Float64}([0.5], [0.5])
```
"""
function convert(::Type{Hyperrectangle}, x::Interval{N, IN}) where {N, IN <: AbstractInterval{N}}
    return Hyperrectangle(low=[low(x)], high=[high(x)])
end

"""
    convert(::Type{HPolytope}, H::AbstractHyperrectangle{N}) where {N}

Converts a hyperrectangular set to a polytope in constraint representation.

### Input

- `HPolytope` -- type used for dispatch
- `H`         -- hyperrectangular set

### Output

A polytope in constraint representation.
"""
function convert(::Type{HPolytope}, H::AbstractHyperrectangle{N}) where {N}
    return HPolytope{N}(constraints_list(H))
end
