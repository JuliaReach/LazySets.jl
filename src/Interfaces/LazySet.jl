import Base: ==, ≈, copy, eltype
import Random.rand

export LazySet,
       basetype,
       ρ, support_function,
       σ, support_vector,
       dim,
       norm,
       radius,
       diameter,
       an_element,
       isbounded,
       neutral,
       absorbing,
       tosimplehrep,
       isuniversal,
       translate,
       affine_map,
       reflect,
       is_interior_point,
       isoperation,
       isoperationtype,
       isequivalent,
       isconvextype,
       area,
       surface,
       singleton_list,
       concretize,
       constraints,
       vertices

"""
    LazySet{N}

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ·x + (1-λ)·y ∈ S``.

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every concrete `LazySet` must define the following functions:
- `σ(d::AbstractVector{N}, S::LazySet{N}) where {N<:Real}` -- the support vector
    of `S` in a given direction `d`; note that the numeric type `N` of `d` and
    `S` must be identical; for some set types `N` may be more restrictive than
    `Real`
- `dim(S::LazySet)` -- the ambient dimension of `S`

The function
- `ρ(d::AbstractVector{N}, S::LazySet{N}) where {N<:Real}` -- the support
    function of `S` in a given direction `d`; note that the numeric type `N` of
    `d` and `S` must be identical; for some set types `N` may be more
    restrictive than `Real`
is optional because there is a fallback implementation relying on `σ`.
However, for unbounded sets (which includes most lazy set types) this fallback
cannot be used and an explicit method must be implemented.

The subtypes of `LazySet` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySet, false)
14-element Array{Any,1}:
 AbstractAffineMap
 AbstractCentrallySymmetric
 AbstractPolyhedron
 Bloating
 CachedMinkowskiSumArray
 CartesianProduct
 CartesianProductArray
 ConvexHull
 ConvexHullArray
 EmptySet
 Intersection
 IntersectionArray
 MinkowskiSum
 MinkowskiSumArray
```

If we only consider *concrete* subtypes, then:

```jldoctest; setup = :(using LazySets: subtypes)
julia> concrete_subtypes = subtypes(LazySet, true);

julia> length(concrete_subtypes)
41

julia> println.(concrete_subtypes);
AffineMap
Ball1
Ball2
BallInf
Ballp
Bloating
CachedMinkowskiSumArray
CartesianProduct
CartesianProductArray
ConvexHull
ConvexHullArray
Ellipsoid
EmptySet
ExponentialMap
ExponentialProjectionMap
HParallelotope
HPolygon
HPolygonOpt
HPolyhedron
HPolytope
HalfSpace
Hyperplane
Hyperrectangle
Intersection
IntersectionArray
Interval
Line2D
LineSegment
LinearMap
Line{N,VN} where VN<:AbstractArray{N,1} where N<:Real
MinkowskiSum
MinkowskiSumArray
ResetMap
Singleton
SymmetricIntervalHull
Translation
Universe
VPolygon
VPolytope
ZeroSet
Zonotope
```
"""
abstract type LazySet{N} end


# --- common LazySet functions ---

"""
    eltype(::Type{<:LazySet{N}}) where {N}

Return the numeric type (`N`) of the given set type.

### Input

- `T` -- set type, used for dispatch

### Output

The numeric type of `T`.
"""
eltype(::Type{<:LazySet{N}}) where {N} = N

"""
    eltype(::LazySet{N}) where {N}

Return the numeric type (`N`) of the given set.

### Input

- `X` -- set instance, used for dispatch

### Output

The numeric type of `X`.
"""
eltype(::LazySet{N}) where {N} = N

"""
    basetype(T::Type{<:LazySet})

Return the base type of the given set type (i.e., without type parameters).

### Input

- `T` -- set type, used for dispatch

### Output

The base type of `T`.
"""
basetype(T::Type{<:LazySet}) = Base.typename(T).wrapper

"""
    basetype(S::LazySet)

Return the base type of the given set (i.e., without type parameters).

### Input

- `S` -- set instance, used for dispatch

### Output

The base type of `S`.

### Examples

```jldoctest
julia> z = rand(Zonotope);

julia> basetype(z)
Zonotope

julia> basetype(z + z)
MinkowskiSum

julia> basetype(LinearMap(rand(2, 2), z + z))
LinearMap
```
"""
basetype(S::LazySet) = Base.typename(typeof(S)).wrapper

"""
    ρ(d::AbstractVector{N}, S::LazySet{N}) where {N<:Real}

Evaluate the support function of a set in a given direction.

### Input

- `d` -- direction
- `S` -- convex set

### Output

The support function of the set `S` for the direction `d`.

### Notes

The numeric type of the direction and the set must be identical.
"""
function ρ(d::AbstractVector{N}, S::LazySet{N}) where {N<:Real}
    return dot(d, σ(d, S))
end

"""
    support_function

Alias for the support function ρ.
"""
const support_function = ρ

"""
    σ

Function to compute the support vector σ.
"""
function σ end

"""
    support_vector

Alias for the support vector σ.
"""
const support_vector = σ

"""
    isbounded(S::LazySet)

Determine whether a set is bounded.

### Input

- `S`         -- set
- `algorithm` -- (optional, default: `"support_function"`) algorithm choice,
                 possible options are `"support_function"` and `"stiemke"`

### Output

`true` iff the set is bounded.

### Algorithm

See the documentation of [`_isbounded_unit_dimensions`](@ref) or
[`_isbounded_stiemke`](@ref) for details.
"""
function isbounded(S::LazySet; algorithm="support_function")
    if algorithm == "support_function"
        return _isbounded_unit_dimensions(S)
    elseif algorithm == "stiemke"
        return _isbounded_stiemke(S)
    else
        throw(ArgumentError("unknown algorithm $algorithm"))
    end
end

"""
    _isbounded_unit_dimensions(S::LazySet{N}) where {N<:Real}

Determine whether a set is bounded in each unit dimension.

### Input

- `S` -- set

### Output

`true` iff the set is bounded in each unit dimension.

### Algorithm

This function performs ``2n`` support function checks, where ``n`` is the
ambient dimension of `S`.
"""
function _isbounded_unit_dimensions(S::LazySet{N}) where {N<:Real}
    n = dim(S)
    @inbounds for i in 1:n
        for o in [one(N), -one(N)]
            d = SingleEntryVector(i, n, o)
            if ρ(d, S) == N(Inf)
                return false
            end
        end
    end
    return true
end

"""
    norm(S::LazySet, [p]::Real=Inf)

Return the norm of a convex set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(S::LazySet, p::Real=Inf)
    if p == Inf
        return norm(Approximations.ballinf_approximation(S), p)
    else
        error("the norm for this value of p=$p is not implemented")
    end
end

"""
    radius(S::LazySet, [p]::Real=Inf)

Return the radius of a convex set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
function radius(S::LazySet, p::Real=Inf)
    if p == Inf
        return radius(Approximations.ballinf_approximation(S), p)
    else
        error("the radius for this value of p=$p is not implemented")
    end
end

"""
    diameter(S::LazySet, [p]::Real=Inf)

Return the diameter of a convex set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(S::LazySet, p::Real=Inf)
    return radius(S, p) * 2
end

"""
    affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector; kwargs...)

Compute a concrete affine map.

### Input

- `M` -- linear map
- `X` -- convex set
- `v` -- translation vector

### Output

A set representing the affine map of `X`.

### Algorithm

The implementation applies the functions `linear_map` and `translate`.
"""
function affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector; kwargs...)
    return translate(linear_map(M, X; kwargs...), v)
end

"""
    an_element(S::LazySet{N}) where {N<:Real}

Return some element of a convex set.

### Input

- `S` -- convex set

### Output

An element of a convex set.

### Algorithm

An element of the set is obtained by evaluating its support vector along
direction ``[1, 0, …, 0]``.
"""
function an_element(S::LazySet{N}) where {N<:Real}
    e₁ = SingleEntryVector(1, dim(S), one(N))
    return σ(e₁, S)
end

"""
    ==(X::LazySet, Y::LazySet)

Return whether two LazySets of the same type are exactly equal.

### Input

- `X` -- any `LazySet`
- `Y` -- another `LazySet` of the same type as `X`

### Output

- `true` iff `X` is equal to `Y`.

### Notes

The check is purely syntactic and the sets need to have the same base type.
For instance, `X::VPolytope == Y::HPolytope` returns `false` even if `X` and `Y`
represent the same polytope.
However `X::HPolytope{Int64} == Y::HPolytope{Float64}` is a valid comparison.

### Algorithm

We recursively compare the fields of `X` and `Y` until a mismatch is found.

### Examples

```jldoctest
julia> HalfSpace([1], 1) == HalfSpace([1], 1)
true

julia> HalfSpace([1], 1) == HalfSpace([1.0], 1.0)
true

julia> Ball1([0.], 1.) == Ball2([0.], 1.)
false
```
"""
function ==(X::LazySet, Y::LazySet)
    # if the common supertype of X and Y is abstract, they cannot be compared
    if isabstracttype(promote_type(typeof(X), typeof(Y)))
        return false
    end

    for f in fieldnames(typeof(X))
        if getfield(X, f) != getfield(Y, f)
            return false
        end
    end

    return true
end
"""
    ≈(X::LazySet, Y::LazySet)

Return whether two LazySets of the same type are approximately equal.

### Input

- `X` -- any `LazySet`
- `Y` -- another `LazySet` of the same type as `X`

### Output

- `true` iff `X` is equal to `Y`.

### Notes

The check is purely syntactic and the sets need to have the same base type.
For instance, `X::VPolytope ≈ Y::HPolytope` returns `false` even if `X` and `Y`
represent the same polytope.
However `X::HPolytope{Int64} ≈ Y::HPolytope{Float64}` is a valid comparison.

### Algorithm

We recursively compare the fields of `X` and `Y` until a mismatch is found.

### Examples

```jldoctest
julia> HalfSpace([1], 1) ≈ HalfSpace([1], 1)
true

julia> HalfSpace([1], 1) ≈ HalfSpace([1.00000001], 0.99999999)
true

julia> HalfSpace([1], 1) ≈ HalfSpace([1.0], 1.0)
true

julia> Ball1([0.], 1.) ≈ Ball2([0.], 1.)
false
```
"""
function ≈(X::LazySet, Y::LazySet)
    # if the common supertype of X and Y is abstract, they cannot be compared
    if isabstracttype(promote_type(typeof(X), typeof(Y)))
        return false
    end

    for f in fieldnames(typeof(X))
        if !_isapprox(getfield(X, f), getfield(Y, f))
            return false
        end
    end

    return true
end

# hook into random API
function rand(rng::AbstractRNG, ::SamplerType{T}) where T<:LazySet
    rand(T, rng=rng)
end

"""
    copy(S::LazySet)

Return a deep copy of the given set by copying its values recursively.

### Input

- `S` -- any `LazySet`

### Output

A copy of `S`.

### Notes

This function performs a `deepcopy` of each field in `S`, resulting in a
completely independent object. See the documentation of `?deepcopy` for further
details.
"""
copy(S::LazySet) = deepcopy(S)

"""
    tosimplehrep(S::LazySet)

Return the simple H-representation ``Ax ≤ b`` of a set from its list of linear
constraints.

### Input

- `S` -- set

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.

### Notes

This function only works for sets that can be represented exactly by a finite
list of linear constraints.
This fallback implementation relies on `constraints_list(S)`.
"""
tosimplehrep(S::LazySet) = tosimplehrep(constraints_list(S))

"""
    reflect(P::LazySet{N}) where {N<:Real}

Concrete reflection of a convex set `P`, resulting in the reflected set `-P`.

### Note

This function requires that the list of constraints of the set `P` is
available, i.e. such that it can be written as
``P = \\{z ∈ ℝⁿ: ⋂ sᵢᵀz ≤ rᵢ, i = 1, ..., N\\}.``

This function can be used to implement the alternative definition of the
Minkowski Difference, which writes as
```math
A ⊖ B = \\{a − b | a ∈ A, b ∈ B\\} = A ⊕ (-B)
```
by calling `minkowski_sum(A, reflect(B))`.
"""
function reflect(P::LazySet)
    @assert applicable(constraints_list, P)  "this function " *
        "requires that the list of constraints is available, but it is not; " *
        "if the set is bounded, try overapproximating with an `HPolytope` first"

    F, g = tosimplehrep(P)
    T = isbounded(P) ? HPolytope : HPolyhedron
    return T(-F, g)
end

"""
    isuniversal(X::LazySet{N}, [witness]::Bool=false) where {N<:Real}

Check whether a given convex set is universal, and otherwise optionally compute
a witness.

### Input

- `X`       -- convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X`` is universal
* If `witness` option is activated:
  * `(true, [])` iff ``X`` is universal
  * `(false, v)` iff ``X`` is not universal and ``v ∉ X``

### Notes

This is a naive fallback implementation.
"""
function isuniversal(X::LazySet{N}, witness::Bool=false) where {N<:Real}
    if isbounded(X)
        result = false
    else
        error("cannot determine universality of the set")
    end

    if result
        return witness ? (true, N[]) : true
    elseif witness
        error("witness production is currently not supported")
    else
        return false
    end
end

"""
    is_interior_point(d::AbstractVector{N}, P::LazySet{N};
                      p=N(Inf), ε=_rtol(N)) where {N<:Real}

Check if the point `d` is contained in the interior of the convex set `P`.

### Input

- `d`  -- point
- `P`  -- set
- `p`  -- (optional; default: `N(Inf)`) norm of the ball used to apply the error
          tolerance
- `ε`  -- (optional; default: `_rtol(N)`) error tolerance of check

### Output

Boolean which indicates if the point `d` is contained in `P`.

### Algorithm

The implementation checks if a `Ballp` of norm `p` with center `d` and radius
`ε` is contained in the set `P`.
This is a numerical check for `d ∈ interior(P)` with error tolerance `ε`.
"""
function is_interior_point(d::AbstractVector{N}, P::LazySet{N};
                           p=N(Inf), ε=_rtol(N)) where {N<:Real}
    return Ballp(p, d, ε) ⊆ P
end

"""
    plot_recipe(X::LazySet{N}, [ε]::N=N(PLOT_PRECISION)) where {N<:Real}

Convert a convex set to a pair `(x, y)` of points for plotting.

### Input

- `X` -- convex set
- `ε` -- (optional, default: `PLOT_PRECISION`) approximation error bound

### Output

A pair `(x, y)` of points that can be plotted.

### Algorithm

We first assert that `X` is bounded.

One-dimensional sets are converted to an `Interval`.
We do not support three-dimensional or higher-dimensional sets at the moment.

For two-dimensional sets, we first compute a polygonal overapproximation.
The second argument, `ε`, corresponds to the error in Hausdorff distance between
the overapproximating set and `X`.
The default value `PLOT_PRECISION` is chosen such that the unit ball in the
2-norm is approximated with reasonable accuracy.
On the other hand, if you only want to produce a fast box-overapproximation of
`X`, pass `ε=Inf`.
Finally, we use the plot recipe for polygons.
"""
function plot_recipe(X::LazySet{N}, ε::N=N(PLOT_PRECISION)) where {N<:Real}
    @assert dim(X) <= 2 "cannot plot a $(dim(X))-dimensional $(typeof(X))"
    @assert isbounded(X) "cannot plot an unbounded $(typeof(X))"

    if dim(X) == 1
        Y = convert(Interval, X)
    else
        Y = overapproximate(X, ε)
    end
    return plot_recipe(Y, ε)
end

"""
    isoperation(X::LazySet)

Check whether the given `LazySet` is an instance of a set operation or not.

### Input

- `X` -- a `LazySet`

### Output

`true` if `X` is an instance of a set-based operation and `false` otherwise.

### Notes

The fallback implementation returns whether the set type of the input is an
operation or not using `isoperationtype`.

See also [`isoperationtype(X::Type{<:LazySet})`](@ref).

### Examples

```jldoctest
julia> B = BallInf([0.0, 0.0], 1.0);

julia> isoperation(B)
false

julia> isoperation(B ⊕ B)
true
```
"""
function isoperation(X::LazySet)
    return isoperationtype(typeof(X))
end

isoperation(::Type{<:LazySet}) = error("`isoperation` cannot be applied to a set " *
                                       "type; use `isoperationtype` instead")

"""
    isoperationtype(X::Type{<:LazySet})

Check whether the given `LazySet` type is an operation or not.

### Input

- `X` -- subtype of `LazySet`

### Output

`true` if the given set type is a set-based operation and `false` otherwise.

### Notes

The fallback for this function returns an error that `isoperationtype` is not
implemented. Subtypes of `LazySet` should dispatch on this function as required.

See also [`isoperation(X<:LazySet)`](@ref).

### Examples

```jldoctest
julia> isoperationtype(BallInf)
false

julia> isoperationtype(LinearMap)
true
```
"""
function isoperationtype(X::Type{<:LazySet})
    error("`isoperationtype` is not implemented for type $X")
end

isoperationtype(::LazySet) = error("`isoperationtype` cannot be applied to " *
                                   "a set instance; use `isoperation` instead")

"""
    isequivalent(X::LazySet, Y::LazySet)

Return whether two LazySets are equal in the mathematical sense, i.e. equivalent.

### Input

- `X` -- any `LazySet`
- `Y` -- another `LazySet`

### Output

`true` iff `X` is equivalent to `Y`.

## Algorithm

First, the check `X == Y` is performed which returns `true` if and only if the given sets are of the same type,
and have the same values (modulo floating-point tolerance). Otherwise, the double inclusion check `X ⊆ Y && Y ⊆ X` is
used.

### Examples

```jldoctest
julia> X = BallInf([0.1, 0.2], 0.3);

julia> Y = convert(HPolytope, X);

julia> X == Y
false

julia> isequivalent(X, Y)
true
```
"""
function isequivalent(X::LazySet, Y::LazySet)
    if X ≈ Y
        return true
    end
    return X ⊆ Y && Y ⊆ X
end

"""
    isconvextype(X::Type{<:LazySet})

Check whether the given `LazySet` type is convex.

### Input

- `X` -- subtype of `LazySet`

### Output

`true` if the given set type is guaranteed to be convex by using only type
information, and `false` otherwise.

### Notes

Since this operation only acts on types (not on values), it can return false
negatives, i.e. there may be instances where the set is convex, even though the
answer of this function is `false`. The examples below illustrate this point.

### Examples

A ball in the infinity norm is always convex, hence we get:

```jldoctest convex_types
julia> isconvextype(BallInf)
true
```

For instance, the union (`UnionSet`) of two sets may in general be either convex
or not, since convexity cannot be decided by just using type information.
Hence, `isconvextype` returns `false` if `X` is `Type{<:UnionSet}`.

```jldoctest convex_types
julia> isconvextype(UnionSet)
false
```

However, the type parameters from the set operations allow to decide convexity
in some cases, by falling back to the convexity of the type of its arguments.
Consider for instance the lazy intersection. The intersection of two convex sets
is always convex, hence we can get:

```jldoctest convex_types
julia> isconvextype(Intersection{Float64, BallInf{Float64}, BallInf{Float64}})
true
```
"""
isconvextype(X::Type{<:LazySet}) = false

"""
    surface(X::LazySet{N}) where {N}

Compute the surface area of a set.

### Input

- `X` -- set

### Output

A number representing the surface area of `X`.
"""
function surface(X::LazySet{N}) where {N}
    if dim(X) == 2
        return area(X)
    else
        throw(ArgumentError("the surface function is only implemented for " *
                    "two-dimensional sets, but the given set is $(dim(X))-dimensional"))
    end
end

"""
    area(X::LazySet{N}) where {N}

Compute the area of a two-dimensional polytopic set using the Shoelace formula.

### Input

- `X` -- two-dimensional set

### Output

A number representing the area of `X`.

### Notes

This algorithm is applicable to any lazy set `X` such that its list of vertices,
`vertices_list`, can be computed.

### Algorithm

Let `m` be the number of vertices of `X`. The following instances are considered:

- `m = 0, 1, 2`: the output is zero.
- `m = 3`: the triangle case is computed using the Shoelace formula with 3 points.
- `m = 4`: the quadrilateral case is obtained by the factored version of the Shoelace
           formula with 4 points.

Otherwise, the general Shoelace formula is used; for detals see the wikipedia
article [Shoelace formula](https://en.wikipedia.org/wiki/Shoelace_formula).
"""
function area(X::LazySet{N}) where {N}
    @assert dim(X) == 2 "this function only applies to two-dimensional sets, " *
    "but the given set is $(dim(X))-dimensional"

    Xpoly = convert(VPolygon, X)
    vlist = vertices_list(Xpoly)
    m = length(vlist)

    if m <= 2
        return zero(N)
    end

    if m == 3 # triangle
        res = _area_triangle(vlist)

    elseif m == 4 # quadrilateral
        res = _area_quadrilateral(vlist)

    else # general case
        res = _area_polygon(vlist)
    end

    return res
end

function _area_triangle(v::Vector{VN}) where {N, VN<:AbstractVector{N}}
    A = v[1]
    B = v[2]
    C = v[3]
    res = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2])
    return abs(res/2)
end

function _area_quadrilateral(v::Vector{VN}) where {N, VN<:AbstractVector{N}}
    A = v[1]
    B = v[2]
    C = v[3]
    D = v[4]
    res = A[1] * (B[2] - D[2]) + B[1] * (C[2] - A[2]) + C[1] * (D[2] - B[2]) + D[1] * (A[2] - C[2])
    return abs(res/2)
end

function _area_polygon(v::Vector{VN}) where {N, VN<:AbstractVector{N}}
    m = length(v)
    @inbounds res = v[m][1] * v[1][2] - v[1][1] * v[m][2]
    for i in 1:m-1
        @inbounds res += v[i][1] * v[i+1][2] - v[i+1][1] * v[i][2]
    end
    return abs(res/2)
end

"""
    singleton_list(P::LazySet)

Return the vertices of a polytopic set as a list of singletons.

### Input

- `P`                 -- polytopic set

### Output

The list of vertices of `P`, as `Singleton`.

### Notes

This function relies on `vertices_list`, which raises an error if the set is
not polytopic (e.g., unbounded).
"""
function singleton_list(P::LazySet)
    return [Singleton(x) for x in vertices_list(P)]
end

"""
    concretize(X::LazySet)

Construct a concrete representation of a (possibly lazy) set.

### Input

- `X` -- set

### Output

A concrete representation of `X` (as far as possible).

### Notes

Since not every lazy set has a concrete set representation in this library, the
result may be partially lazy.
"""
function concretize(X::LazySet)
    return X
end

"""
    constraints(X::LazySet)

Construct an iterator over the constraints of a polyhedral set.

### Input

- `X` -- polyhedral set

### Output

An iterator over the constraints of `X`.
"""
function constraints(X::LazySet)
    return _constraints_fallback(X)
end

"""
    vertices(X::LazySet)

Construct an iterator over the vertices of a polyhedral set.

### Input

- `X` -- polyhedral set

### Output

An iterator over the vertices of `X`.
"""
function vertices(X::LazySet)
    return _vertices_fallback(X)
end

function _constraints_fallback(X::LazySet)
    return VectorIterator(constraints_list(X))
end

function _vertices_fallback(X::LazySet)
    return VectorIterator(vertices_list(X))
end

# =========================
# Code requiring MiniQhull
# =========================

function load_delaunay_MiniQhull()
return quote

import .MiniQhull: delaunay
export delaunay

"""
    delaunay(X::LazySet)

Compute the Delaunay triangulation of the given convex set.

### Input

- `X` -- set

### Output

A union of polytopes in vertex representation.

### Notes

This function requires that you have properly installed the package
[MiniQhull.jl](https://github.com/gridap/MiniQhull.jl), including the library
[Qhull](http://www.qhull.org/).

The method works in arbitrary dimension and the requirement is that the list of
vertices of `X` can be obtained.
"""
function delaunay(X::LazySet)
    n = dim(X)
    v = vertices_list(X)
    m = length(v)
    coordinates = vcat(v...)
    connectivity_matrix = delaunay(n, m, coordinates)
    nelements = size(connectivity_matrix, 2)
    elements = [VPolytope(v[connectivity_matrix[:, i]]) for i in 1:nelements]
    return UnionSetArray(elements)
end

end end  # load_delaunay_MiniQhull
