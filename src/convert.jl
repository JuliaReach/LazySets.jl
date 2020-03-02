import Base.convert

#= conversion between set types =#

# convert methods for identity (no-ops)
for T in subtypes(LazySet, true)
    @eval begin
        Base.convert(::Type{$T}, X::$T) = X
    end
end

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

"""
    convert(::Type{HPOLYGON}, P::HPolytope{N}) where
        {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from 2D polytope in H-representation to polygon in H-representation.

### Input

- `type`  -- target type
- `P`     -- source polytope (must be 2D)
- `prune` -- (optional, default: `true`) flag for removing redundant
             constraints in the end
### Output

The 2D polytope represented as polygon.
"""
function convert(::Type{HPOLYGON}, P::HPolytope{N};
                 prune::Bool=true) where {N<:Real, HPOLYGON<:AbstractHPolygon}
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    H = HPOLYGON{N}()
    for ci in constraints_list(P)
        addconstraint!(H, ci; prune=prune)
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
    convert(::Type{Zonotope}, Z::AbstractZonotope)

Converts a zonotopic set to a zonotope.

### Input

- `Zonotope`
- `H` -- zonotopic set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, Z::AbstractZonotope)
    return Zonotope(center(Z), genmat(Z))
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}) where {N<:Real,
        HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}

Converts the cartesian product of two hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `cp`       -- cartesian product of two hyperrectangular sets

### Output

This method falls back to the conversion of the cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}
                ) where {N<:Real, HN1<:AbstractHyperrectangle{N},
                         HN2<:AbstractHyperrectangle{N}}
    return convert(Zonotope, convert(Hyperrectangle, cp))
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN})
        where {N<:Real, HN<:AbstractHyperrectangle{N}}

Converts the cartesian product array of hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `cpa`      -- cartesian product array of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method falls back to the conversion of the cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN}
                ) where {N<:Real, HN<:AbstractHyperrectangle{N}}
    return convert(Zonotope, convert(Hyperrectangle, cpa))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, ZN}
           ) where {N, ZN<:AbstractZonotope{N}}

Converts the lazy linear map of a zonotopic set to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of a zonotopic set

### Output

A zonotope.

### Algorithm

This method first applies the (concrete) linear map to the zonotopic set and
then converts the result to a `Zonotope` type.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, ZN}
                ) where {N, ZN<:AbstractZonotope{N}}
    return convert(Zonotope, linear_map(S.M, S.X))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
           ) where {N, HN1<:AbstractHyperrectangle{N},
                    HN2<:AbstractHyperrectangle{N}}

Converts the lazy linear map of the cartesian product of two hyperrectangular
sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of the cartesian product of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the cartesian product to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
                ) where {N, HN1<:AbstractHyperrectangle{N},
                         HN2<:AbstractHyperrectangle{N}}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope},S::LinearMap{N, CartesianProductArray{N, HN}}
           ) where {N, HN<:AbstractHyperrectangle{N}}

Converts the lazy linear map of the cartesian product of a finite number of
hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of a `CartesianProductArray` of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the cartesian product array to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, CartesianProductArray{N, HN}}
                ) where {N, HN<:AbstractHyperrectangle{N}}
    return linear_map(S.M, convert(Zonotope, S.X))
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
    constraints_list = Vector{LinearConstraint{N, Vector{N}}}(undef, 3)
    o = one(N)
    z = zero(N)
    v = element(S)
    constraints_list[1] = LinearConstraint([o, o], v[1] + v[2])
    constraints_list[2] = LinearConstraint([-o, z], -v[1])
    constraints_list[3] = LinearConstraint([z, -o], -v[2])
    return HPOLYGON(constraints_list)
end

"""
    convert(::Type{HPOLYGON}, L::LineSegment{N}
          ) where {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from line segment to polygon in H-representation.

### Input

- `type`  -- target type
- `L`     -- line segment
- `prune` -- (optional, default: `false`) flag for removing redundant
             constraints in the end
### Output

A flat polygon in constraint representation with the minimal number of
constraints (four).
"""
function convert(::Type{HPOLYGON}, L::LineSegment{N}) where {N<:Real, HPOLYGON<:AbstractHPolygon}
    H = HPOLYGON{N}()
    c = halfspace_left(L.p, L.q)
    addconstraint!(H, c; prune=false)
    addconstraint!(H, LinearConstraint(-c.a, -c.b); prune=false)
    line_dir = L.q - L.p
    c = LinearConstraint(line_dir, dot(L.q, line_dir))
    addconstraint!(H, c; prune=false)
    line_dir = -line_dir
    addconstraint!(H, LinearConstraint(line_dir, dot(L.p, line_dir)); prune=false)
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
Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}([0.5], [0.5])
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
julia> convert(Interval, Hyperrectangle([0.5], [0.5]))
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```
"""
function convert(::Type{Interval}, H::AbstractHyperrectangle)
    @assert dim(H) == 1 "cannot convert a $(dim(H))-dimensional $(typeof(H)) to `Interval`"
    return Interval(low(H)[1], high(H)[1])
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
    @assert dim(S) == 1 "cannot convert a $(dim(H))-dimensional $(typeof(S)) to `Interval`"
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
            cp::CartesianProduct{N, HN1, HN2}) where {N<:Real, HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}

Converts the cartesian product of two hyperrectangular sets to a single hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

The result is obtained by concatenating the center and radius of each hyperrectangle.
This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cp::CartesianProduct{N, HN1, HN2}) where {N<:Real, HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}
     X, Y = cp.X, cp.Y
     c = vcat(center(X), center(Y))
     r = vcat(radius_hyperrectangle(X), radius_hyperrectangle(Y))
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
    convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N<:Real}

Converts a two-dimensional hyperrectangle to the cartesian product of two
intervals.

### Input

- `CartesianProduct` -- type used for dispatch
- `H`                -- hyperrectangle

### Output

The cartesian product of two intervals.
"""
function convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N<:Real}
    @assert dim(H) == 2 "the hyperrectangle must be two-dimensional to convert it to " *
            "the cartesian product of two intervals, but it is $(dim(H))-dimensional; " *
            "consider converting it to a 'CartesianProductArray{$N, Interval{$N}}' instead"
    Ix = Interval(low(H, 1), high(H, 1))
    Iy =  Interval(low(H, 2), high(H, 2))
    return CartesianProduct(Ix, Iy)
end

"""
    convert(::Type{CartesianProductArray{N, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N<:Real}

Converts a hyperrectangle to the cartesian product array of intervals.

### Input

- `CartesianProductArray` -- type used for dispatch
- `H`                     -- hyperrectangle

### Output

The cartesian product of a finite number of intervals.
"""
function convert(::Type{CartesianProductArray{N, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N<:Real}
    Iarray = [Interval(low(H, i), high(H, i)) for i in 1:dim(H)]
    return CartesianProductArray(Iarray)
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

The cartesian product is obtained by:

- Concatenating the centers of each input zonotope.
- Arranging the generators in block-diagional fashion, and filled with zeros
  in the off-diagonal; for this reason, the generator matrix of the returned
  zonotope is built as a sparse matrix.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, <:Zonotope{N}, <:Zonotope{N}}) where {N<:Real}
    Z1, Z2 = cp.X, cp.Y
    c = vcat(Z1.center, Z2.center)
    G = blockdiag(sparse(Z1.generators), sparse(Z2.generators))
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
    convert(::Type{IntervalArithmetic.Interval}, x::Interval)

Converts a `LazySets` interval to an `Interval` from `IntervalArithmetic`.

### Input

- `Interval` -- type used for dispatch, from `IntervalArithmetic`
- `x`        -- interval (`LazySets.Interval`)

### Output

An `IntervalArithmetic.Interval`.
"""
function convert(::Type{IntervalArithmetic.Interval}, x::Interval)
    return IntervalArithmetic.interval(min(x), max(x))
end

"""
    convert(::Type{Interval}, x::IntervalArithmetic.Interval)

Converts an `Interval` from `IntervalArithmetic` to an interval in `LazySets`.

### Input

- `Interval` -- type used for dispatch
- `x`        -- interval (`IntervalArithmetic.Interval`)

### Output

A `LazySets.Interval`.
"""
function convert(::Type{Interval}, x::IntervalArithmetic.Interval)
    return Interval(IntervalArithmetic.inf(x), IntervalArithmetic.sup(x))
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

### Notes

`IntervalArithmetic.IntervalBox` uses *static* vectors to store each component
interval, hence the resulting `Hyperrectangle` has its center and radius represented
as a static vector (`SArray`).
"""
function convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)
    low_IB = IntervalArithmetic.inf.(IB)
    high_IB = IntervalArithmetic.sup.(IB)
    return Hyperrectangle(low=low_IB, high=high_IB)
end

"""
    convert(::Type{Hyperrectangle}, r::Rectification{N, AH})
        where {N<:Real, AH<:AbstractHyperrectangle{N}}

Converts a rectification of a hyperrectangle to a hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `r`              -- rectification of a hyperrectangle

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 r::Rectification{N, AH}) where {N<:Real,
                                                 AH<:AbstractHyperrectangle{N}}
    return rectify(r.X)
end

"""
    convert(::Type{Interval},
            r::Rectification{N, IN}) where {N<:Real, IN<:Interval{N}}

Converts a rectification of an interval to an interval.

### Input

- `Interval` -- type used for dispatch
- `r`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval},
                 r::Rectification{N, IN}) where {N<:Real, IN<:Interval{N}}
    return Interval(rectify([min(r.X), max(r.X)]))
end

"""
    convert(::Type{VPolytope},
            X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}

Converts the convex hull array of singletons to a polytope in V-representation.

### Input

- `VPolytope` -- type used for dispatch
- `X`         -- convex hull array of singletons

### Output

A polytope in vertex representation.
"""
function convert(::Type{VPolytope},
                 X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    return VPolytope(vertices_list(X))
end

"""
    convert(::Type{VPolygon},
            X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}

Converts the convex hull array of singletons to a polygon in V-representation.

### Input

- `VPolygon`  -- type used for dispatch
- `X`         -- convex hull array of singletons

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon},
                 X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    @assert dim(X) == 2
    return VPolygon(vertices_list(X))
end

"""
    convert(::Type{MinkowskiSumArray},
            X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}

Converts the Minkowski sum of a Minkowski sum array to a Minkowski sum array.

### Input

- `MinkowskiSumArray`  -- type used for dispatch
- `X`                  -- Minkowski sum of a Minkowski sum array

### Output

A Minkowski sum array.
"""
function convert(::Type{MinkowskiSumArray},
                 X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}
    return MinkowskiSumArray(vcat(X.X, X.Y.array))
end

"""
    convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval{N}}

Converts the Minkowski sum of two intervals into an interval.

### Input

- `Interval` -- type used for dispatch
- `x`        -- Minkowski sum of a pair of intervals

### Output

An interval.
"""
function convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval{N}}
    return minkowski_sum(x.X, x.Y)
end
