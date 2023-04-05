import Base.convert

# convert methods for identity (no-ops)
for T in subtypes(LazySet, true)
    @eval begin
        Base.convert(::Type{$T}, X::$T) = X
    end
end

"""
    convert(T::Type{HPOLYGON}, P::VPolygon) where {HPOLYGON<:AbstractHPolygon}

Convert a polygon in vertex representation to a polygon in constraint
representation.

### Input

- `HPOLYGON` -- target type
- `P`        -- polygon in vertex representation

### Output

A polygon in constraint representation.
"""
function convert(T::Type{HPOLYGON},
                 P::VPolygon) where {HPOLYGON<:AbstractHPolygon}
    return tohrep(P, T)
end

"""
    convert(::Type{HPOLYGON}, X::LazySet; [check_boundedness]::Bool=true,
            prune::Bool=true) where {HPOLYGON<:AbstractHPolygon}

Convert a two-dimensional polytope to a polygon in vertex representation.

### Input

- `HPOLYGON`          -- target type
- `X`                 -- two-dimensional polytope
- `check_boundedness` -- (optional, default `!isboundedtype(typeof(X))`) if
                         `true` check whether the set `X` is bounded before
                         creating the polygon
- `prune`             -- (optional, default: `true`) flag for removing redundant
                         constraints in the end

### Output

A polygon in constraint representation.

### Algorithm

We compute the list of constraints of `X`, then instantiate the polygon.
"""
function convert(::Type{HPOLYGON}, X::LazySet;
                 check_boundedness::Bool=!isboundedtype(typeof(X)),
                 prune::Bool=true) where {HPOLYGON<:AbstractHPolygon}
    @assert dim(X) == 2 "set must be two-dimensional for conversion, but it " *
        "has dimension $(dim(X))"
    PT = basetype(HPOLYGON)
    if check_boundedness && !isbounded(X)
        throw(ArgumentError("expected a bounded set for conversion to `$PT`"))
    end
    return PT(constraints_list(X); prune=prune)
end

"""
    convert(::Type{VPolygon}, P::AbstractHPolygon)

Convert a polygon in constraint representation to a polygon in vertex
representation.

### Input

- `VPolygon` -- target type
- `P`        -- polygon in constraint representation

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon}, P::AbstractHPolygon)
    return tovrep(P)
end

"""
    convert(::Type{VPolygon}, X::LazySet)

Convert a two-dimensional polytopic set to a polygon in vertex representation.

### Input

- `VPolygon` -- target type
- `X`        -- two-dimensional polytopic set

### Output

The 2D set represented as a polygon.

### Algorithm

This method uses `vertices_list`.
"""
function convert(::Type{VPolygon}, X::LazySet)
    @assert dim(X) == 2 "set must be two-dimensional for conversion"
    return VPolygon(vertices_list(X))
end

"""
    convert(::Type{HPolytope}, X::LazySet)

Convert a polytopic set to a polytope in constraint representation.

### Input

- `HPolytope` -- target type
- `X`         -- polytopic set

### Output

The given polytope represented as a polytope in constraint representation.

### Algorithm

This method uses `constraints_list`.
"""
function convert(::Type{HPolytope}, X::LazySet)
    if !isboundedtype(typeof(X)) || !is_polyhedral(X)
        error("conversion to `HPolytope` requires a polytopic set")
    end
    return HPolytope(constraints_list(X))
end

"""
    convert(::Type{HPolyhedron}, X::LazySet)

Convert a polyhedral set to a polyhedron in constraint representation.

### Input

- `HPolyhedron` -- target type
- `X`           -- polyhedral set

### Output

The given set represented as a polyhedron in constraint representation.

### Algorithm

This method uses `constraints_list`.
"""
function convert(::Type{HPolyhedron}, X::LazySet)
    if !is_polyhedral(X)
        error("conversion to `HPolyhedron` requires a polyhedral set")
    end
    return HPolyhedron(constraints_list(X))
end

# conversion of a lazyset to a polyhedron with dense normal vectors
function convert(::Type{HPolyhedron{N, Vector{N}}}, X::LazySet) where {N}
    if !is_polyhedral(X)
        error("conversion to `HPolyhedron` requires a polyhedral set")
    end
    return HPolyhedron([HalfSpace(Vector(c.a), c.b) for c in constraints(X)])
end

"""
    convert(::Type{VPolytope}, X::LazySet; [prune]::Bool=true)

Convert a polytopic set to a polytope in vertex representation.

### Input

- `VPolytope` -- target type
- `X`         -- polytopic set
- `prune`     -- (optional, default: `true`) option to remove redundant vertices

### Output

The given set represented as a polytope in vertex representation.

### Algorithm

This method uses `vertices_list`. Use the option `prune` to select whether to
remove redundant vertices before constructing the polytope.
"""
function convert(::Type{VPolytope}, X::LazySet; prune::Bool=true)
    return VPolytope(vertices_list(X, prune=prune))
end

# fast conversion from a 2D hyperrectangular set to a zonotope
function _convert_2D(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    c = center(H)
    rx = radius_hyperrectangle(H, 1)
    ry = radius_hyperrectangle(H, 2)
    G = _genmat_2D(c, rx, ry)
    return Zonotope(c, G)
end

@inline function _genmat_2D(c::AbstractVector{N}, rx, ry) where {N}
    flat_x = isapproxzero(rx)
    flat_y = isapproxzero(ry)
    ncols = !flat_x + !flat_y
    G = Matrix{N}(undef, 2, ncols)
    if !flat_x
        @inbounds begin G[1] = rx; G[2] = zero(N) end
        if !flat_y
            @inbounds begin G[3] = zero(N); G[4] = ry end
        end
    elseif !flat_y
        @inbounds begin G[1] = zero(N); G[2] = ry end
    end
    return G
end

function load_genmat_2D_static()
return quote

@inline function _genmat_2D(c::SVector{L, N}, rx, ry) where {L, N}
    flat_x = isapproxzero(rx)
    flat_y = isapproxzero(ry)
    if !flat_x && !flat_y
        G = SMatrix{2, 2, N, 4}(rx, zero(N), zero(N), ry)
    elseif !flat_x && flat_y
        G = SMatrix{2, 1, N, 2}(rx, zero(N))
    elseif flat_x && !flat_y
        G = SMatrix{2, 1, N, 2}(zero(N), ry)
    else
        G = SMatrix{2, 0, N, 0}()
    end
    return G
end

# this function is type-stable but doesn't prune the generators according
# to flat dimensions of H
function _convert_2D_static(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    c = center(H)
    rx = radius_hyperrectangle(H, 1)
    ry = radius_hyperrectangle(H, 2)
    G = SMatrix{2, 2, N, 4}(rx, zero(N), zero(N), ry)
    return Zonotope(c, G)
end

function _convert_static(::Type{Zonotope}, H::Hyperrectangle{N, <:SVector, <:SVector}) where {N}
    return Zonotope(center(H), _genmat_static(H))
end

end end  # quote / load_genmat_2D_static

"""
    convert(::Type{Zonotope}, Z::AbstractZonotope)

Convert a zonotopic set to a zonotope.

### Input

- `Zonotope` -- target type
- `H`        -- zonotopic set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, Z::AbstractZonotope)
    return _convert_zonotope_fallback(Z)
end

function convert(::Type{Zonotope}, H::AbstractHyperrectangle)
    dim(H) == 2 && return _convert_2D(Zonotope, H)
    return _convert_zonotope_fallback(H)
end

_convert_zonotope_fallback(Z) = Zonotope(center(Z), genmat(Z))

function convert(::Type{Singleton}, cp::CartesianProduct{N, S1, S2}
                ) where {N, S1<:AbstractSingleton, S2<:AbstractSingleton}
    return Singleton(vcat(element(cp.X), element(cp.Y)))
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}) where {N,
            HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}

Convert the Cartesian product of two hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cp`       -- Cartesian product of two hyperrectangular sets

### Output

This method falls back to the conversion of the Cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}
                ) where {N, HN1<:AbstractHyperrectangle,
                            HN2<:AbstractHyperrectangle}
    return convert(Zonotope, convert(Hyperrectangle, cp))
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN})
        where {N, HN<:AbstractHyperrectangle}

Convert the Cartesian product array of hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cpa`      -- Cartesian product array of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method falls back to the conversion of the Cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN}
                ) where {N, HN<:AbstractHyperrectangle}
    return convert(Zonotope, convert(Hyperrectangle, cpa))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, ZN})
        where {N, ZN<:AbstractZonotope}

Convert the lazy linear map of a zonotopic set to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a zonotopic set

### Output

A zonotope.

### Algorithm

This method first applies the (concrete) linear map to the zonotopic set and
then converts the result to a `Zonotope` type.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, ZN}
                ) where {N, ZN<:AbstractZonotope}
    return convert(Zonotope, linear_map(S.M, S.X))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
           ) where {N, HN1<:AbstractHyperrectangle,
                    HN2<:AbstractHyperrectangle}

Convert the lazy linear map of the Cartesian product of two hyperrectangular
sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of the Cartesian product of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
                ) where {N, HN1<:AbstractHyperrectangle,
                         HN2<:AbstractHyperrectangle}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope},S::LinearMap{N, CartesianProductArray{N, HN}})
        where {N, HN<:AbstractHyperrectangle}

Convert the lazy linear map of the Cartesian product of a finite number of
hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- linear map of a `CartesianProductArray` of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the Cartesian product array to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, CartesianProductArray{N, HN}}
                ) where {N, HN<:AbstractHyperrectangle}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
           ) where {N, HPOLYGON<:AbstractHPolygon}

Convert a singleton to a polygon in constraint representation.

### Input

- `HPOLYGON` -- target type
- `S`        -- singleton

### Output

A polygon in constraint representation with the minimal number of constraints
(three).
"""
function convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
                ) where {N, HPOLYGON<:AbstractHPolygon}
    constraints_list = Vector{HalfSpace{N, Vector{N}}}(undef, 3)
    o = one(N)
    z = zero(N)
    v = element(S)
    constraints_list[1] = HalfSpace([o, o], v[1] + v[2])
    constraints_list[2] = HalfSpace([-o, z], -v[1])
    constraints_list[3] = HalfSpace([z, -o], -v[2])
    return HPOLYGON(constraints_list)
end

"""
    convert(::Type{HPOLYGON}, L::LineSegment{N}
          ) where {N, HPOLYGON<:AbstractHPolygon}

Convert a line segment to a polygon in constraint representation.

### Input

- `HPOLYGON` -- target type
- `L`        -- line segment
- `prune`    -- (optional, default: `false`) flag for removing redundant
                constraints in the end
### Output

A flat polygon in constraint representation with the minimal number of
constraints (four).
"""
function convert(::Type{HPOLYGON}, L::LineSegment{N}) where {N, HPOLYGON<:AbstractHPolygon}
    H = HPOLYGON{N}()
    c = halfspace_left(L.p, L.q)
    addconstraint!(H, c; prune=false)
    addconstraint!(H, HalfSpace(-c.a, -c.b); prune=false)
    line_dir = L.q - L.p
    c = HalfSpace(line_dir, dot(L.q, line_dir))
    addconstraint!(H, c; prune=false)
    line_dir = -line_dir
    addconstraint!(H, HalfSpace(line_dir, dot(L.p, line_dir)); prune=false)
    return H
end

"""
    convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)

Convert a hyperrectangular set to a hyperrectangle.

### Input

- `Hyperrectangle` -- hyperrectangle target type
- `H`              -- hyperrectangular set

### Output

A hyperrectangle.

### Examples

```jldoctest
julia> convert(Hyperrectangle, Interval(0.0, 1.0))
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([0.5], [0.5])
```
"""
function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

"""
    convert(::Type{Interval}, X::LazySet)

Convert a one-dimensional convex set to an interval.

### Input

- `Interval` -- interval target type
- `X`        -- one-dimensional convex set

### Output

An interval.

### Examples

```jldoctest
julia> convert(Interval, Hyperrectangle([0.5], [0.5]))
Interval{Float64, IntervalArithmetic.Interval{Float64}}([0, 1])
```
"""
function convert(::Type{Interval}, X::LazySet)
    @assert dim(X) == 1 "cannot convert a $(dim(X))-dimensional set to an " *
                        "`Interval`"
    if !isconvextype(typeof(X))
        error("this implementation requires a convex set")
    end

    l, h = extrema(X, 1)
    return Interval(l, h)
end

"""
    convert(::Type{Hyperrectangle}, cpa::CartesianProductArray{N, HN})
        where {N, HN<:AbstractHyperrectangle}

Convert the Cartesian product of a finite number of hyperrectangular sets to
a single hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle}, cpa::CartesianProductArray{N, HN}
                ) where {N, HN<:AbstractHyperrectangle}
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
    convert(::Type{Hyperrectangle}, cp::CartesianProduct{N, HN1, HN2})
        where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}

Convert the Cartesian product of two hyperrectangular sets to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

The result is obtained by concatenating the center and radius of each
hyperrectangle. This implementation uses the `center` and
`radius_hyperrectangle` methods.
"""
function convert(::Type{Hyperrectangle}, cp::CartesianProduct{N, HN1, HN2}
                ) where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}
     X, Y = cp.X, cp.Y
     c = vcat(center(X), center(Y))
     r = vcat(radius_hyperrectangle(X), radius_hyperrectangle(Y))
     return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, IN}) where {N, IN<:Interval}

Convert the Cartesian product of a finite number of intervals to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of intervals

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `min` and `max` methods of `Interval` to reduce
the allocations and improve performance (see LazySets#1143).
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N, IN}) where {N, IN<:Interval}
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
            H::AbstractHyperrectangle{N}) where {N}

Convert a two-dimensional hyperrectangle to the Cartesian product of two
intervals.

### Input

- `CartesianProduct` -- target type
- `H`                -- hyperrectangle

### Output

The Cartesian product of two intervals.
"""
function convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    @assert dim(H) == 2 "the hyperrectangle must be two-dimensional to " *
            "convert it to the Cartesian product of two intervals, but it is " *
            "$(dim(H))-dimensional; consider converting it to a " *
            "`CartesianProductArray{$N, Interval{$N}}` instead"
    Ix = Interval(low(H, 1), high(H, 1))
    Iy =  Interval(low(H, 2), high(H, 2))
    return CartesianProduct(Ix, Iy)
end

"""
    convert(::Type{CartesianProductArray{N, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N}

Convert a hyperrectangle to the Cartesian product array of intervals.

### Input

- `CartesianProductArray` -- target type
- `H`                     -- hyperrectangle

### Output

The Cartesian product of a finite number of intervals.
"""
function convert(::Type{CartesianProductArray{N, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    Iarray = [Interval(low(H, i), high(H, i)) for i in 1:dim(H)]
    return CartesianProductArray(Iarray)
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, ZN1, ZN2}
           ) where {N, ZN1<:AbstractZonotope, ZN2<:AbstractZonotope}

Convert the Cartesian product of two zonotopic sets to a new zonotope.

### Input

- `Zonotope` -- target type
- `S`        -- Cartesian product of two zonotopic sets

### Output

A zonotope.

### Algorithm

The Cartesian product is obtained by:

- Concatenating the centers of each input zonotope.
- Arranging the generators in block-diagonal fashion, and filled with zeros in
  the off-diagonal; for this reason, the generator matrix of the returned
  zonotope is built as a sparse matrix.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, ZN1, ZN2}
                ) where {N, ZN1<:AbstractZonotope, ZN2<:AbstractZonotope}
    Z1, Z2 = cp.X, cp.Y
    c = vcat(center(Z1), center(Z2))
    G = blockdiag(sparse(genmat(Z1)), sparse(genmat(Z2)))
    return Zonotope(c, G)
end

"""
    convert(::Type{IntervalArithmetic.Interval}, X::LazySet)

Convert a convex set to an `Interval` from `IntervalArithmetic`.

### Input

- `Interval` -- target type, from `IntervalArithmetic`
- `X`        -- convex set

### Output

An `IntervalArithmetic.Interval`.
"""
function convert(::Type{IA.Interval}, X::LazySet)
    @assert dim(X) == 1 "cannot convert a $(dim(X))-dimensional set to an " *
                        "interval"
    return convert(Interval, X).dat
end

"""
    convert(::Type{Interval}, x::IntervalArithmetic.Interval)

Convert an `Interval` from `IntervalArithmetic` to an `Interval` in `LazySets`.

### Input

- `Interval` -- target type
- `x`        -- interval (`IntervalArithmetic.Interval`)

### Output

A `LazySets.Interval`.
"""
function convert(::Type{Interval}, x::IA.Interval)
    return Interval(x)
end

"""
    convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)

Convert a hyperrectangular set to an `IntervalBox` from `IntervalArithmetic`.

### Input

- `IntervalBox` -- target type
- `H`           -- hyperrectangular set

### Output

An `IntervalBox`.
"""
function convert(::Type{IA.IntervalBox}, H::AbstractHyperrectangle)
    return IA.IntervalBox(IA.interval.(low(H), high(H)))
end

"""
    convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)

Convert an `IntervalBox` from `IntervalArithmetic` to a hyperrectangular set.

### Input

- `Hyperrectangle` -- target type
- `IB`             -- interval box

### Output

A `Hyperrectangle`.

### Notes

`IntervalArithmetic.IntervalBox` uses *static* vectors to store each component
interval; hence the resulting `Hyperrectangle` has its center and radius
represented as a static vector (`SArray`).
"""
function convert(::Type{Hyperrectangle}, IB::IA.IntervalBox)
    low_IB = IA.inf.(IB)
    high_IB = IA.sup.(IB)
    return Hyperrectangle(low=low_IB, high=high_IB)
end

# method for Interval
function convert(::Type{Hyperrectangle}, I::IA.Interval)
    low_I = [IA.inf(I)]
    high_I = [IA.sup(I)]
    return Hyperrectangle(low=low_I, high=high_I)
end

"""
    convert(::Type{Hyperrectangle}, r::Rectification{N, AH})
        where {N, AH<:AbstractHyperrectangle}

Convert a rectification of a hyperrectangle to a hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `r`              -- rectification of a hyperrectangle

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle}, r::Rectification{N, AH}
                ) where {N, AH<:AbstractHyperrectangle}
    return rectify(r.X)
end

"""
    convert(::Type{Interval}, r::Rectification{N, IN}) where {N, IN<:Interval}

Convert a rectification of an interval to an interval.

### Input

- `Interval` -- target type
- `r`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval},
                 r::Rectification{N, IN}) where {N, IN<:Interval}
    return Interval(rectify([min(r.X), max(r.X)]))
end

"""
    convert(::Type{MinkowskiSumArray},
            X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}

Convert the Minkowski sum of a Minkowski sum array to a Minkowski sum array.

### Input

- `MinkowskiSumArray`  -- target type
- `X`                  -- Minkowski sum of a Minkowski sum array

### Output

A Minkowski sum array.
"""
function convert(::Type{MinkowskiSumArray},
                 X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}
    return MinkowskiSumArray(vcat(X.X, X.Y.array))
end

"""
    convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval}

Convert the Minkowski sum of two intervals to an interval.

### Input

- `Interval` -- target type
- `x`        -- Minkowski sum of two intervals

### Output

An interval.
"""
function convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval}
    return minkowski_sum(x.X, x.Y)
end

# convert to concrete Vector representation
function convert(::Type{HalfSpace{N, Vector{N}}},
                 hs::HalfSpace{N, <:AbstractVector{N}}) where {N}
    return HalfSpace(Vector(hs.a), hs.b)
end

function convert(::Type{Zonotope},
                 am::AbstractAffineMap{N, <:AbstractZonotope{N}}) where {N}
    Z1 = convert(Zonotope, linear_map(matrix(am), set(am)))
    translate!(Z1, vector(am))
    return Z1
end

"""
    convert(::Type{HParallelotope}, Z::AbstractZonotope{N}) where {N}

Convert a zonotopic set of order one to a parallelotope in constraint
representation.

### Input

- `HParallelotope` -- target type
- `Z`              -- zonotopic set of order one

### Output

A parallelotope in constraint representation.

### Notes

This function requires that the list of constraints of `Z` are obtained in
the particular order returned from the `constraints_list` function of a
`Zonotope`. Hence it first converts `Z` to a `Zonotope`.
"""
function convert(::Type{HParallelotope}, Z::AbstractZonotope{N}) where {N}
    @assert order(Z) == 1 "cannot convert a zonotope that is not of order 1 " *
                          "to a parallelotope"
    n = dim(Z)

    constraints = _constraints_list_zonotope(Z)

    D = Matrix{N}(undef, n, n)
    c = Vector{N}(undef, 2n)
    j = 1
    @inbounds for i in 1:n
        D[i, :] = constraints[j].a
        c[i] = constraints[j].b
        c[i+n] = constraints[j+1].b
        j += 2
    end
    return HParallelotope(D, c; check_consistency=false)
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, AZ})
        where {N, AZ<:AbstractZonotope}

Convert a Cartesian product array of zonotopic sets to a zonotope.

### Input

- `Zonotope` -- target type
- `cpa`       -- Cartesian product array of zonotopic sets

### Output

A zonotope with sparse matrix representation.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N, AZ}
                ) where {N, AZ<:AbstractZonotope}
    arr = array(cpa)
    c = reduce(vcat, center.(arr))
    G = reduce(blockdiag, sparse.(genmat.(arr)))
    return Zonotope(c, G)
end

# make a copy of the constraints
convert(::Type{HPolytope}, P::HPolyhedron) = HPolytope(copy(constraints_list(P)))
convert(::Type{HPolyhedron}, P::HPolytope) = HPolyhedron(copy(constraints_list(P)))

for T in [HPolygon, HPolygonOpt, HPolytope, HPolyhedron]
@eval begin
    function convert(::Type{$T}, P::Intersection)
        clist = vcat(constraints_list(P.X), constraints_list(P.Y))
        return ($T)(clist)
    end

    function convert(::Type{$T}, P::IntersectionArray)
        clist = reduce(vcat, constraints_list.(array(P)))
        return ($T)(clist)
    end
end
end

"""
    convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}

Convert a polyhedral set to a star set represented as a lazy affine map.

### Input

- `STAR` -- target type
- `P`    -- polyhedral set

### Output

A star set.
"""
function convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N)*I, n, n)
    return AffineMap(V, P, c)
end

"""
    convert(::Type{STAR}, X::Star)

Convert a star set to its equivalent representation as a lazy affine map.

### Input

- `STAR` -- target type
- `X`    -- star set

### Output

A star set.
"""
function convert(::Type{STAR}, X::Star)
    return AffineMap(X.V, X.P, X.c)
end

"""
    convert(::Type{Star}, P::AbstractPolyhedron{N}) where {N}

Convert a polyhedral set to a star set.

### Input

- `Star` -- target type
- `P`    -- polyhedral set

### Output

A star set.
"""
function convert(::Type{Star}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N)*I, n, n)
    return Star(c, V, P)
end

function convert(::Type{Hyperplane}, P::HPolyhedron; skip_check::Bool=false)
    # check that the number of constraints is fine
    if !skip_check && !is_hyperplanar(P)
        throw(ArgumentError("the polyhedron is not hyperplanar: $P"))
    end

    # construct hyperplane from first constraint
    c1 = @inbounds first(constraints_list(P))
    return Hyperplane(c1.a, c1.b)
end

"""
    convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)

Convert a zonotope to a simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `Z`                              -- zonotopic set

### Output

A simple sparse polynomial zonotope.

### Algorithm

This method implements Proposition 3 in [1].

[1] Kochdumper, Althoff. *Sparse polynomial zonotopes - a novel set
representation for reachability analysis*. 2021
"""
function convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    n = ngens(Z)
    E = Matrix(1 * I, n, n)

    return SimpleSparsePolynomialZonotope(c, G, E)
end

"""
    convert(::Type{SimpleSparsePolynomialZonotope},
            SPZ::SparsePolynomialZonotope)

Convert a sparse polynomial zonotope to simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `SPZ`                            -- sparse polynomial zonotope

### Output

A simple sparse polynomial zonotope.
"""
function convert(::Type{SimpleSparsePolynomialZonotope},
                 SPZ::SparsePolynomialZonotope)
    c = center(SPZ)
    G = hcat(genmat_dep(SPZ), genmat_indep(SPZ))
    n = ngens_indep(SPZ)
    E = cat(expmat(SPZ), Matrix(1 * I, n, n), dims=(1, 2))

    return SimpleSparsePolynomialZonotope(c, G, E)
end

"""
    convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}

Convert a zonotope to sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `Z`                        -- zonotopic set

### Output

A sparse polynomial zonotope.
"""
function convert(::Type{SparsePolynomialZonotope},
                 Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    G = genmat(Z)
    n = ngens(Z)
    E = Matrix(1 * I, n, n)
    idx = uniqueID(n)
    GI = zeros(N, dim(Z), 0)

    return SparsePolynomialZonotope(c, G, GI, E, idx)
end

"""
    convert(::Type{SparsePolynomialZonotope},
            SSPZ::SimpleSparsePolynomialZonotope)

Convert a simple sparse polynomial zonotope to a sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `SSPZ`                     -- simple sparse polynomial zonotope

### Output

A sparse polynomial zonotope.
"""
function convert(::Type{SparsePolynomialZonotope},
                 SSPZ::SimpleSparsePolynomialZonotope{N}) where {N}
    c = center(SSPZ)
    G = genmat(SSPZ)
    E = expmat(SSPZ)
    n = ngens(SSPZ)
    idx = uniqueID(n)
    GI = zeros(N, dim(SSPZ), 0)

    return SparsePolynomialZonotope(c, G, GI, E, idx)
end
