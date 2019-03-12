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

First the list of constraints of `P` is computed, then the corresponding
`HPolytope` is created.
"""
function convert(::Type{HPolytope}, P::AbstractPolytope)
    return HPolytope(constraints_list(P))
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
        addconstraint!(H, ci; prune=false)
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
    addconstraint!(H, c; prune=false)
    addconstraint!(H, LinearConstraint(-c.a, -c.b); prune=false)
    line_dir = L.q - L.p
    c = LinearConstraint(line_dir, dot(L.q, line_dir))
    addconstraint!(H, c; prune=false)
    line_dir = -line_dir
    addconstraint!(H, LinearConstraint(line_dir, dot(L.p, line_dir));
                   prune=false)
    return H
end

"""
    convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)

Convert a hyperrectangular set to a hyperrectangle.

### Input

- `Hyperrectangle` -- hyperrectangle type, used for dispatch
- `H`              -- hyperrectangular set

### Output

A hyperrectangle.

### Examples

```jldoctest
julia> convert(Hyperrectangle, Interval(0.0, 1.0))
Hyperrectangle{Float64}([0.5], [0.5])
```
"""
function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

"""
    convert(::Type{Interval}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to an interval.

### Input

- `Interval` -- interval type, used for dispatch
- `H`        -- hyperrectangular set

### Output

An interval.

### Examples

```jldoctest
julia> convert(Interval, Hyperrectangle{Float64}([0.5], [0.5]))
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```
"""
function convert(::Type{Interval}, H::AbstractHyperrectangle)
    @assert dim(H) == 1 "can only convert a one-dimensional $(typeof(H)) to `Interval`"
    return Interval([low(H); high(H)])
end

"""
    convert(::Type{Interval}, S::LazySet{N}) where {N<:Real}

Converts a convex set to an interval.

### Input

- `Interval` -- interval type, used for dispatch
- `S`        -- convex set

### Output

An interval.
"""
function convert(::Type{Interval}, S::LazySet{N}) where {N<:Real}
    @assert dim(S) == 1 "can only convert a one-dimensional $(typeof(S)) to `Interval`"
    return Interval(-ρ(N[-1], S), ρ(N[1], S))
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, HN}) where {N<:Real, HN<:AbstractHyperrectangle{N}}

Converts the cartesian product of a finite number of hyperrectangular sets to
a single hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product array of hyperrectangular set 

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N, HN}) where {N<:Real, HN<:AbstractHyperrectangle{N}}
     n = dim(cpa)
     c = Vector{N}(undef, n)
     r = Vector{N}(undef, n)
     i = 1
     @inbounds for block_set in array(cpa)
         j = i + dim(block_set) - 1
         c[i:j] = center(block_set)
         r[i:j] = radius_hyperrectangle(block_set)
         i = j + 1
     end
     return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, Interval{N}}) where {N<:Real}

Converts the cartesian product of a finite number of intervals to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product array of intervals 

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `min` and `max` methods of `Interval` to reduce
the allocatons and improve performance (see LazySets#1143).
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N, Interval{N}}) where {N<:Real}
     # since the sets are intervals, the dimension of cpa is its length
     n = length(array(cpa))
     l = Vector{N}(undef, n)
     h = Vector{N}(undef, n)
     @inbounds for (i, Ii) in enumerate(array(cpa))
         l[i] = min(Ii)
         h[i] = max(Ii)
     end
     return Hyperrectangle(low=l, high=h)
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, Zonotope{N}, Zonotope{N}}) where {N<:Real}

Converts the cartesian product of two zonotopes to a new zonotope.

### Input

- `Zonotope` -- type used for dispatch
- `S`        -- cartesian product of two zonotopes

### Output

A zonotope.

### Algorithm

This implementation concatenates the centers of each input zonotope.
The resulting generator matrix is such that the generators for each element of the
cartesian product are added along the diagonal.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, Zonotope{N}, Zonotope{N}}) where {N<:Real}
    Z1, Z2 = cp.X, cp.Y
    n1, p1 = size(Z1.generators)
    n2, p2 = size(Z2.generators)
    c = vcat(Z1.center, Z2.center)
    G = [Z1.generators zeros(n1, p2);
         zeros(n2, p1) Z2.generators]
    return Zonotope(c, G)
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

"""
    convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to an `IntervalBox` from `IntervalArithmetic`.

### Input

- `IntervalBox` -- type used for dispatch
- `H`           -- hyperrectangular set

### Output

An `IntervalBox`.
"""
function convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)
    return IntervalArithmetic.IntervalBox(IntervalArithmetic.interval.(low(H), high(H)))
end

"""
    convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)

Converts an `IntervalBox` from `IntervalArithmetic` to a hyperrectangular set.

### Input

- `Hyperrectangle` -- type used for dispatch
- `IB`             -- interval box

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)
    low_IB = Vector(IntervalArithmetic.inf.(IB))    # TODO: temprary conversion, see #1214
    high_IB = Vector(IntervalArithmetic.sup.(IB))   # TODO: temprary conversion, see #1214
    return Hyperrectangle(low=low_IB, high=high_IB)
end
