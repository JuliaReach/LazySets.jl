import Base.convert

#= conversion between set types =#

"""
    convert(::Type{HPOLYGON1},
            P::HPOLYGON2) where {HPOLYGON1<:AbstractHPolygon,
                                 HPOLYGON2<:AbstractHPolygon}

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
                 P::HPOLYGON2) where {HPOLYGON1<:AbstractHPolygon,
                                      HPOLYGON2<:AbstractHPolygon}
    if P isa HPOLYGON1
        return P
    end
    return HPOLYGON1(P.constraints)
end

"""
    convert(T::Type{HPOLYGON}, P::VPolygon) where {HPOLYGON<:AbstractHPolygon}

Converts a polygon in vertex representation to a polygon in constraint representation.

### Input

- `HPOLYGON` -- type used for dispatch
- `P`        -- polygon in vertex representation

### Output

A polygon in constraint representation.
"""
function convert(T::Type{HPOLYGON}, P::VPolygon) where {HPOLYGON<:AbstractHPolygon}
    return tohrep(P, T)
end

"""
    convert(::Type{VPolygon}, P::AbstractHPolygon)

Converts a polygon in constraint representation to a polygon in vertex representation.

### Input

- `VPolygon` -- type used for dispatch
- `P`        -- polygon in constraint representation

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon}, P::AbstractHPolygon)
    return tovrep(P)
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
    convert(::Type{VPolytope}, P::AbstractPolytope)

Convert polytopic type to polytope in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The set `P` represented as a `VPolytope`.
"""
function convert(::Type{VPolytope}, P::AbstractPolytope)
    return VPolytope(vertices_list(P))
end

"""
    convert(::Type{VPolygon}, P::AbstractPolytope)

Convert polytopic set to polygon in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The 2D polytope represented as a polygon.
"""
function convert(::Type{VPolygon}, P::AbstractPolytope)
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    return VPolygon(vertices_list(P))
end

"""
    convert(::Type{HPolytope}, P::AbstractPolytope)

Convert a polytopic set to a polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The given polytope represented as a polytope in constraint representation.

### Algorithm

``P`` is first converted to a polytope in V-representation.
Then, the conversion method to a polytope in H-representation is invoked.
This conversion may require the `Polyhedra` library.
"""
function convert(::Type{HPolytope}, P::AbstractPolytope)
    return convert(HPolytope, convert(VPolytope, P))
end

"""
    convert(::Type{HPolyhedron}, P::AbstractPolytope)

Convert a polytopic set to a polyhedron in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The given polytope represented as a polyhedron in constraint representation.
"""
function convert(::Type{HPolyhedron}, P::AbstractPolytope)
    return HPolyhedron(constraints_list(P))
end

"""
    convert(::Type{HPolytope}, P::VPolytope)

Convert from polytope in V-representation to polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The polytope in the dual representation.

### Algorithm

The `tohrep` function is invoked. It requires the `Polyhedra` package.
"""
function convert(::Type{HPolytope}, P::VPolytope)
    return tohrep(P)
end

"""
    convert(::Type{VPolytope}, P::HPolytope)

Convert from polytope in H-representation to polytope in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The polytope in the dual representation.

### Algorithm

The `tovrep` function is invoked. It requires the `Polyhedra` package.
"""
function convert(::Type{VPolytope}, P::HPolytope)
    return tovrep(P)
end

# no-op's
convert(::Type{HPolytope}, P::HPolytope) = P
convert(::Type{VPolytope}, P::VPolytope) = P

"""
    convert(::Type{HPOLYGON}, P::HPolytope{N}) where
        {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from 2D polytope in H-representation to polygon in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope (must be 2D)

### Output

The 2D polytope represented as polygon.
"""
function convert(::Type{HPOLYGON},
                 P::HPolytope{N}) where {N<:Real, HPOLYGON<:AbstractHPolygon}
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    H = HPOLYGON{N}()
    for ci in constraints_list(P)
        # guarantee that the edges are correctly sorted for storage
        addconstraint!(H, ci)
    end
    return H
end

"""
    convert(::Type{Zonotope}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to a zonotope.

### Input

- `Zonotope`
- `H` -- hyperrectangular set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, H::AbstractHyperrectangle)
    return Zonotope(center(H), Diagonal(radius_hyperrectangle(H)))
end

"""
    convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
           ) where {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from singleton to polygon in H-representation.

### Input

- `type` -- target type
- `S`    -- singleton

### Output

A polygon in constraint representation with the minimal number of constraints
(three).
"""
function convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
                ) where {N<:Real, HPOLYGON<:AbstractHPolygon}
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
          ) where {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from line segment to polygon in H-representation.

### Input

- `type` -- target type
- `L`    -- line segment

### Output

A flat polygon in constraint representation with the minimal number of
constraints (four).
"""
function convert(::Type{HPOLYGON}, L::LineSegment{N}
                ) where {N<:Real, HPOLYGON<:AbstractHPolygon}
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
    convert(::Type{Hyperrectangle}, x::Interval)

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
function convert(::Type{Hyperrectangle}, x::Interval)
    return Hyperrectangle(low=[min(x)], high=[max(x)])
end

"""
    convert(::Type{HPolytope}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to a polytope in constraint representation.

### Input

- `HPolytope` -- type used for dispatch
- `H`         -- hyperrectangular set

### Output

A polytope in constraint representation.
"""
function convert(::Type{HPolytope},
                 H::AbstractHyperrectangle)
    return HPolytope(constraints_list(H))
end

"""
    convert(::Type{HPOLYGON}, H::AbstractHyperrectangle) where
        {HPOLYGON<:AbstractHPolygon}

Converts a hyperrectangular set to a polygon in constraint representation.

### Input

- `HPOLYGON`  -- type used for dispatch
- `H`         -- hyperrectangular set

### Output

A polygon in constraint representation.
"""
function convert(X::Type{HPOLYGON},
                 H::AbstractHyperrectangle) where {HPOLYGON<:AbstractHPolygon}
    @assert dim(H) == 2 "cannot convert a $(dim(H))-dimensional " *
        "hyperrectangle into a two-dimensional polygon"
    return HPOLYGON(constraints_list(H))
end
