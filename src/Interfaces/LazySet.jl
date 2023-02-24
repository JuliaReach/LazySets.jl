import Base: extrema, ==, ≈, copy, eltype, rationalize
import Random.rand

export LazySet,
       low,
       high,
       basetype,
       ρ, support_function,
       σ, support_vector,
       complement,
       dim,
       norm,
       radius,
       diameter,
       an_element,
       isbounded,
       isboundedtype,
       neutral,
       absorbing,
       tosimplehrep,
       isuniversal,
       translate,
       affine_map,
       exponential_map,
       reflect,
       is_interior_point,
       isoperation,
       isoperationtype,
       isequivalent,
       isconvextype,
       is_polyhedral,
       area,
       surface,
       singleton_list,
       concretize,
       constraints,
       vertices,
       project,
       rectify,
       permute,
       chebyshev_center_radius

"""
    LazySet{N}

Abstract type for the set types in LazySets.

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every concrete `LazySet` must define the following method:

- `dim(S::LazySet)` -- the ambient dimension of `S`

While not strictly required, it is useful to define the following method:

- `σ(d::AbstractVector, S::LazySet)` -- the support vector of `S` in a given
    direction `d`

The method

- `ρ(d::AbstractVector, S::LazySet)` -- the support function of `S` in a given
    direction `d`

is optional because there is a fallback implementation relying on `σ`.
However, for potentially unbounded sets (which includes most lazy set types)
this fallback cannot be used and an explicit method must be implemented.

The subtypes of `LazySet` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySet, false)
18-element Vector{Any}:
 AbstractAffineMap
 AbstractPolynomialZonotope
 Bloating
 CachedMinkowskiSumArray
 CartesianProduct
 CartesianProductArray
 Complement
 ConvexSet
 Intersection
 IntersectionArray
 LazySets.AbstractStar
 MinkowskiSum
 MinkowskiSumArray
 Polygon
 QuadraticMap
 Rectification
 UnionSet
 UnionSetArray
```

If we only consider *concrete* subtypes, then:

```jldoctest; setup = :(using LazySets: subtypes)
julia> concrete_subtypes = subtypes(LazySet, true);

julia> length(concrete_subtypes)
54

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
Complement
ConvexHull
ConvexHullArray
DensePolynomialZonotope
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
InverseLinearMap
LazySets.AbstractStar
Line
Line2D
LineSegment
LinearMap
MinkowskiSum
MinkowskiSumArray
Polygon
QuadraticMap
Rectification
ResetMap
RotatedHyperrectangle
SimpleSparsePolynomialZonotope
Singleton
SparsePolynomialZonotope
Star
SymmetricIntervalHull
Translation
UnionSet
UnionSetArray
Universe
VPolygon
VPolytope
ZeroSet
Zonotope

```
"""
abstract type LazySet{N} end

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
negatives, i.e., there may be instances where the set is convex, even though the
answer of this function is `false`. The examples below illustrate this point.

### Examples

A ball in the infinity norm is always convex, hence we get:

```jldoctest convex_types
julia> isconvextype(BallInf)
true
```

For instance, the union (`UnionSet`) of two sets may in general be either convex
or not. Since convexity cannot be decided by just using type information,
`isconvextype` returns `false`.

```jldoctest convex_types
julia> isconvextype(UnionSet)
false
```

However, the type parameters of set operations allow to decide convexity in some
cases by falling back to the convexity information of the type of its arguments.
Consider for instance the lazy intersection. The intersection of two convex sets
is always convex, hence we get:

```jldoctest convex_types
julia> isconvextype(Intersection{Float64, BallInf{Float64}, BallInf{Float64}})
true
```
"""
isconvextype(X::Type{<:LazySet}) = false

# Note: this method cannot be documented due to a bug in Julia
function low(X::LazySet, i::Int)
    return _low(X, i)
end

function _low(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, -one(N))
    return -ρ(d, X)
end

"""
    low(X::LazySet)

Return a vector with the lowest coordinates of the set in each canonical
direction.

### Input

- `X` -- set

### Output

A vector with the lower coordinate of the set in each dimension.

### Notes

See also `low(X::LazySet, i::Int)`.

The result is the lowermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function low(X::LazySet)
    n = dim(X)
    return [low(X, i) for i in 1:n]
end

# Note: this method cannot be documented due to a bug in Julia
function high(X::LazySet, i::Int)
    return _high(X, i)
end

function _high(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, one(N))
    return ρ(d, X)
end

"""
    high(X::LazySet)

Return a vector with the highest coordinate of the set in each canonical
direction.

### Input

- `X` -- set

### Output

A vector with the highest coordinate of the set in each dimension.

### Notes

See also `high(X::LazySet, i::Int)`.

The result is the uppermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function high(X::LazySet)
    n = dim(X)
    return [high(X, i) for i in 1:n]
end

"""
    extrema(X::LazySet, i::Int)

Return the lower and higher coordinate of a set in a given dimension.

### Input

- `X` -- set
- `i` -- dimension of interest

### Output

The lower and higher coordinate of the set in the given dimension.

### Notes

The result is equivalent to `(low(X, i), high(X, i))`, but sometimes it can be
computed more efficiently.

### Algorithm

The bounds are computed with `low` and `high`.
"""
function extrema(X::LazySet, i::Int)
    l = low(X, i)
    h = high(X, i)
    return (l, h)
end

"""
    extrema(X::LazySet)

Return two vectors with the lowest and highest coordinate of `X` in each
canonical direction.

### Input

- `X` -- set

### Output

Two vectors with the lowest and highest coordinates of `X` in each dimension.

### Notes

See also [`extrema(X::LazySet, i::Int)`](@ref).

The result is equivalent to `(low(X), high(X))`, but sometimes it can be
computed more efficiently.

The resulting points are the lowermost and uppermost corners of the box
approximation, so they are not necessarily contained in `X`.

### Algorithm

The bounds are computed with `low` and `high` by default.
"""
function extrema(X::LazySet)
    l = low(X)
    h = high(X)
    return (l, h)
end

"""
    plot_vlist(X::S, ε::Real) where {S<:LazySet}

Return a list of vertices used for plotting a two-dimensional set.

### Input

- `X` -- two-dimensional set
- `ε` -- precision parameter

### Output

A list of vertices of a polygon `P`.
For convex `X`, `P` usually satisfies that the Hausdorff distance to `X` is less
than `ε`.
"""
function plot_vlist(X::S, ε::Real) where {S<:LazySet}
    @assert isconvextype(S) "can only plot convex sets"

    P = overapproximate(X, ε)
    return convex_hull(vertices_list(P))
end

"""
    convex_hull(X::LazySet; kwargs...)

Compute the convex hull of a polytopic set.

### Input

- `X` -- polytopic set

### Output

The set `X` itself if its type indicates that it is convex, or a new set with
the list of the vertices describing the convex hull.

### Algorithm

For non-convex sets, this method relies on the `vertices_list` method.
"""
function convex_hull(X::LazySet; kwargs...)
    if isconvextype(typeof(X))
        return X
    end

    vlist = convex_hull(vertices_list(X); kwargs...)
    return _convex_hull_set(vlist; n=dim(X))
end

"""
    eltype(::Type{<:LazySet{N}}) where {N}

Return the numeric type (`N`) of the given set type.

### Input

- `T` -- set type

### Output

The numeric type of `T`.
"""
eltype(::Type{<:LazySet{N}}) where {N} = N

"""
    eltype(::LazySet{N}) where {N}

Return the numeric type (`N`) of the given set.

### Input

- `X` -- set

### Output

The numeric type of `X`.
"""
eltype(::LazySet{N}) where {N} = N

"""
    basetype(T::Type{<:LazySet})

Return the base type of the given set type (i.e., without type parameters).

### Input

- `T` -- set type

### Output

The base type of `T`.
"""
basetype(T::Type{<:LazySet}) = Base.typename(T).wrapper

"""
    basetype(S::LazySet)

Return the base type of the given set (i.e., without type parameters).

### Input

- `S` -- set

### Output

The base type of `S`.

### Examples

```jldoctest
julia> Z = rand(Zonotope);

julia> basetype(Z)
Zonotope

julia> basetype(Z + Z)
MinkowskiSum

julia> basetype(LinearMap(rand(2, 2), Z + Z))
LinearMap
```
"""
basetype(S::LazySet) = basetype(typeof(S))

"""
    ρ(d::AbstractVector, S::LazySet)

Evaluate the support function of a set in a given direction.

### Input

- `d` -- direction
- `S` -- set

### Output

The evaluation of the support function of the set `S` for the direction `d`.
"""
function ρ(d::AbstractVector, S::LazySet)
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
    isboundedtype(T::Type{<:LazySet})

Check whether a set type only represents bounded sets.

### Input

- `T` -- set type

### Output

`true` if the set type only represents bounded sets.
Note that some sets may still represent an unbounded set even though their type
actually does not (example: [`HPolytope`](@ref), because the construction with
non-bounding linear constraints is allowed).

### Notes

By default this function returns `false`.
All set types that can determine boundedness should override this behavior.
"""
function isboundedtype(T::Type{<:LazySet})
    return false
end

"""
    isbounded(S::LazySet)

Check whether a set is bounded.

### Input

- `S`         -- set
- `algorithm` -- (optional, default: `"support_function"`) algorithm choice,
                 possible options are `"support_function"` and `"stiemke"`

### Output

`true` iff the set is bounded.

### Algorithm

See the documentation of `_isbounded_unit_dimensions` or `_isbounded_stiemke`
for details.
"""
function isbounded(S::LazySet; algorithm="support_function")
    if algorithm == "support_function"
        return _isbounded_unit_dimensions(S)
    elseif algorithm == "stiemke"
        return _isbounded_stiemke(constraints_list(S))
    else
        throw(ArgumentError("unknown algorithm $algorithm"))
    end
end

"""
    _isbounded_unit_dimensions(S::LazySet)

Check whether a set is bounded in each unit dimension.

### Input

- `S` -- set

### Output

`true` iff the set is bounded in each unit dimension.

### Algorithm

This function asks for upper and lower bounds in each ambient dimension.
"""
function _isbounded_unit_dimensions(S::LazySet)
    @inbounds for i in 1:dim(S)
        if isinf(low(S, i)) || isinf(high(S, i))
            return false
        end
    end
    return true
end

"""
    is_polyhedral(S::LazySet)

Trait for polyhedral sets.

### Input

- `S` -- set

### Output

`true` only if the set behaves like an [`AbstractPolyhedron`](@ref).

### Notes

The answer is conservative, i.e., may sometimes be `false` even if the set is
polyhedral.
"""
function is_polyhedral(S::LazySet)
    return false
end

"""
    norm(S::LazySet, [p]::Real=Inf)

Return the norm of a set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(S::LazySet, p::Real=Inf)
    if p == Inf
        return norm(box_approximation(S), p)
    elseif is_polyhedral(S) && isboundedtype(typeof(S))
        return maximum(norm(v, p) for v in vertices_list(S))
    else
        error("the norm for this value of p=$p is not implemented")
    end
end

"""
    radius(S::LazySet, [p]::Real=Inf)

Return the radius of a set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- set
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

Return the diameter of a set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(S::LazySet, p::Real=Inf)
    return radius(S, p) * 2
end

"""
    affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector; kwargs...)

Compute the concrete affine map ``M·X + v``.

### Input

- `M` -- linear map
- `X` -- set
- `v` -- translation vector

### Output

A set representing the affine map ``M·X + v``.

### Algorithm

The implementation applies the functions `linear_map` and `translate`.
"""
function affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector; kwargs...)
    return translate(linear_map(M, X; kwargs...), v)
end

"""
    exponential_map(M::AbstractMatrix, X::LazySet)

Compute the concrete exponential map of `M` and `X`, i.e., `exp(M) * X`.

### Input

- `M` -- matrix
- `X` -- set

### Output

A set representing the exponential map of `M` and `X`.

### Algorithm

The implementation applies the functions `exp` and `linear_map`.
"""
function exponential_map(M::AbstractMatrix, X::LazySet)
    return linear_map(exp(M), X)
end

"""
    an_element(S::LazySet)

Return some element of a set.

### Input

- `S` -- set

### Output

An element of a set.

### Algorithm

An element of the set is obtained by evaluating its support vector along
direction ``[1, 0, …, 0]``.
This may fail for unbounded sets.
"""
function an_element(S::LazySet)
    return _an_element_lazySet(S)
end

function _an_element_lazySet(S::LazySet)
    N = eltype(S)
    e₁ = SingleEntryVector(1, dim(S), one(N))
    return σ(e₁, S)
end

"""
    ==(X::LazySet, Y::LazySet)

Check whether two sets use exactly the same set representation.

### Input

- `X` -- set
- `Y` -- set

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

julia> Ball1([0.0], 1.0) == Ball2([0.0], 1.0)
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

Check whether two sets of the same type are approximately equal.

### Input

- `X` -- set
- `Y` -- set of the same base type as `X`

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

julia> Ball1([0.0], 1.0) ≈ Ball2([0.0], 1.0)
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

Return a copy of a set by copying its values recursively.

### Input

- `S` -- set

### Output

A copy of `S`.

### Notes

This function computes a `copy` of each field in `S`.
See the documentation of `?copy` for further details.
"""
function copy(S::T) where {T<:LazySet}
    args = [copy(getfield(S, f)) for f in fieldnames(T)]
    BT = basetype(S)
    return BT(args...)
end

"""
    tosimplehrep(S::LazySet)

Return the simple constraint representation ``Ax ≤ b`` of a polyhedral set from
its list of linear constraints.

### Input

- `S` -- polyhedral set

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.

### Algorithm

This fallback implementation relies on `constraints_list(S)`.
"""
tosimplehrep(S::LazySet) = tosimplehrep(constraints_list(S))

"""
    reflect(P::LazySet)

Concrete reflection of a set `P`, resulting in the reflected set `-P`.

### Algorithm

This function requires that the list of constraints of the set `P` is
available, i.e., that it can be written as
``P = \\{z ∈ ℝⁿ: ⋂ sᵢᵀz ≤ rᵢ, i = 1, ..., N\\}.``

This function can be used to implement the alternative definition of the
Minkowski Difference
```math
A ⊖ B = \\{a − b | a ∈ A, b ∈ B\\} = A ⊕ (-B)
```
by calling `minkowski_sum(A, reflect(B))`.
"""
function reflect(P::LazySet)
    if !is_polyhedral(P)
        error("this implementation requires a polyhedral set; try " *
              "overapproximating with an `HPolyhedron` first")
    end

    F, g = tosimplehrep(P)
    T = isbounded(P) ? HPolytope : HPolyhedron
    return T(-F, g)
end

"""
    is_interior_point(d::AbstractVector{N}, X::LazySet{N};
                      p=N(Inf), ε=_rtol(N)) where {N}

Check whether the point `d` is contained in the interior of the set `X`.

### Input

- `d`  -- point
- `X`  -- set
- `p`  -- (optional; default: `N(Inf)`) norm of the ball used to apply the error
          tolerance
- `ε`  -- (optional; default: `_rtol(N)`) error tolerance of check

### Output

Boolean which indicates if the point `d` is contained in `X`.

### Algorithm

The implementation checks if a `Ballp` of norm `p` with center `d` and radius
`ε` is contained in the set `X`.
This is a numerical check for `d ∈ interior(X)` with error tolerance `ε`.
"""
function is_interior_point(d::AbstractVector{N}, X::LazySet{N};
                           p=N(Inf), ε=_rtol(N)) where {N}
    return Ballp(p, d, ε) ⊆ X
end

"""
    plot_recipe(X::LazySet, [ε])

Convert a compact convex set to a pair `(x, y)` of points for plotting.

### Input

- `X` -- compact convex set
- `ε` -- approximation-error bound

### Output

A pair `(x, y)` of points that can be plotted.

### Notes

We do not support three-dimensional or higher-dimensional sets at the moment.

### Algorithm

One-dimensional sets are converted to an `Interval`.

For two-dimensional sets, we first compute a polygonal overapproximation.
The second argument, `ε`, corresponds to the error in Hausdorff distance between
the overapproximating set and `X`.
On the other hand, if you only want to produce a fast box-overapproximation of
`X`, pass `ε=Inf`.

Finally, we use the plot recipe for the constructed set (interval or polygon).
"""
function plot_recipe(X::LazySet, ε)
    @assert dim(X) <= 2 "cannot plot a $(dim(X))-dimensional $(typeof(X))"
    @assert isboundedtype(typeof(X)) || isbounded(X) "cannot plot an " *
        "unbounded $(typeof(X))"
    @assert isconvextype(typeof(X)) "can only plot convex sets"

    if dim(X) == 1
        Y = convert(Interval, X)
    else
        Y = overapproximate(X, ε)
    end
    return plot_recipe(Y, ε)
end

"""
    isoperation(X::LazySet)

Check whether a set is an instance of a set operation or not.

### Input

- `X` -- set

### Output

`true` if `X` is an instance of a set-based operation and `false` otherwise.

### Notes

This fallback implementation checks whether the set type of the input is an
operation type using [`isoperationtype(::Type{<:LazySet})`](@ref).

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

# common error
isoperation(::Type{<:LazySet}) = error("`isoperation` cannot be applied to a " *
                                      "set type; use `isoperationtype` instead")

"""
    isoperationtype(X::Type{<:LazySet})

Check whether the given set type is an operation or not.

### Input

- `X` -- set type

### Output

`true` if the given set type is a set operation and `false` otherwise.

### Notes

This fallback implementation returns an error that `isoperationtype` is not
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

# common error
isoperationtype(::LazySet) = error("`isoperationtype` cannot be applied to a " *
                                   "set instance; use `isoperation` instead")

"""
    isequivalent(X::LazySet, Y::LazySet)

Check whether two sets are equal in the mathematical sense, i.e., equivalent.

### Input

- `X` -- set
- `Y` -- set

### Output

`true` iff `X` is equivalent to `Y` (up to some precision).

## Algorithm

First we check `X ≈ Y`, which returns `true` if and only if `X` and `Y` have the
same type and approximately the same values (checked with `LazySets._isapprox`).
If that fails, we check the double inclusion `X ⊆ Y && Y ⊆ X`.

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
    try  # TODO temporary try-catch construct until ≈ is fixed for all set types
        if X ≈ Y
            return true
        end
    catch e
    end
    return _isequivalent_inclusion(X, Y)
end

function _isequivalent_inclusion(X::LazySet, Y::LazySet)
    return X ⊆ Y && Y ⊆ X
end

"""
    surface(X::LazySet)

Compute the surface area of a set.

### Input

- `X` -- set

### Output

A real number representing the surface area of `X`.
"""
function surface(X::LazySet)
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

- `X` -- two-dimensional polytopic set

### Output

A number representing the area of `X`.

### Notes

This algorithm is applicable to any polytopic set `X` whose list of vertices can
be computed via `vertices_list`.

### Algorithm

Let `m` be the number of vertices of `X`. We consider the following instances:

- `m = 0, 1, 2`: the output is zero.
- `m = 3`: the triangle case is solved using the Shoelace formula with 3 points.
- `m = 4`: the quadrilateral case is solved by the factored version of the
           Shoelace formula with 4 points.

Otherwise, the general Shoelace formula is used; for details see the
[Wikipedia page](https://en.wikipedia.org/wiki/Shoelace_formula).
"""
function area(X::LazySet{N}) where {N}
    @assert isconvextype(typeof(X)) "this function requires a convex set"
    @assert dim(X) == 2 "this function only applies to two-dimensional sets, " *
        "but the given set is $(dim(X))-dimensional"

    Xpoly = convert(VPolygon, X)  # sorts vertices
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
    res = A[1] * (B[2] - D[2]) + B[1] * (C[2] - A[2]) + C[1] * (D[2] - B[2]) +
          D[1] * (A[2] - C[2])
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

- `P` -- polytopic set

### Output

A list of the vertices of `P` as `Singleton`s.

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

Construct an iterator over the vertices of a polytopic set.

### Input

- `X` -- polytopic set

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

function load_delaunay_MiniQhull()
return quote

import .MiniQhull: delaunay
export delaunay

"""
    delaunay(X::LazySet)

Compute the Delaunay triangulation of the given polytopic set.

### Input

- `X` -- polytopic set

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

"""
    complement(X::LazySet)

Return the complement of a polyhedral set.

### Input

- `X` -- polyhedral set

### Output

A `UnionSetArray` of half-spaces, i.e., the output is the union of the linear
constraints which are obtained by complementing each constraint of `X`.

### Algorithm

The principle used in this implementation is that for any pair of sets
``(X, Y)`` we have that ``(X ∩ Y)^C = X^C ∪ Y^C``. In particular, we can apply
this rule for each constraint that defines a polyhedral set. Hence the concrete
complement can be represented as the set union of the complement of each
constraint.
"""
function complement(X::LazySet)
    return UnionSetArray(constraints_list(Complement(X)))
end

"""
    project(S::LazySet, block::AbstractVector{Int}, [::Nothing=nothing],
            [n]::Int=dim(S); [kwargs...])

Project a set to a given block by using a concrete linear map.

### Input

- `S`       -- set
- `block`   -- block structure - a vector with the dimensions of interest
- `nothing` -- (default: `nothing`)
- `n`       -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the projection of the set `S` to block `block`.

### Algorithm

We apply the function `linear_map`.
"""
@inline function project(S::LazySet, block::AbstractVector{Int},
                         ::Nothing=nothing, n::Int=dim(S); kwargs...)
    return _project_linear_map(S, block, n; kwargs...)
end

@inline function _project_linear_map(S::LazySet{N}, block::AbstractVector{Int},
                                     n::Int=dim(S); kwargs...) where {N}
    M = projection_matrix(block, n, N)
    return linear_map(M, S)
end

"""
    project(S::LazySet, block::AbstractVector{Int}, set_type::Type{TS},
            [n]::Int=dim(S); [kwargs...]) where {TS<:LazySet}

Project a set to a given block and set type, possibly involving an
overapproximation.

### Input

- `S`        -- set
- `block`    -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n`        -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set of type `set_type` representing an overapproximation of the projection of
`S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set using `overapproximate` and `set_type`.
"""
@inline function project(S::LazySet, block::AbstractVector{Int},
                         set_type::Type{TS}, n::Int=dim(S);
                         kwargs...) where {TS<:LazySet}
    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type)
end

"""
    project(S::LazySet, block::AbstractVector{Int},
            set_type_and_precision::Pair{T, N}, [n]::Int=dim(S);
            [kwargs...]) where {T<:UnionAll, N<:Real}

Project a set to a given block and set type with a certified error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type_and_precision` -- pair `(T, ε)` of a target set type `T` and an
                              error bound `ε` for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Notes

Currently we only support `HPolygon` as set type, which implies that the set
must be two-dimensional.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set with the given error bound `ε`.
"""
@inline function project(S::LazySet, block::AbstractVector{Int},
                         set_type_and_precision::Pair{T, N}, n::Int=dim(S);
                         kwargs...) where {T<:UnionAll, N<:Real}
    set_type, ε = set_type_and_precision
    @assert length(block) == 2 && set_type == HPolygon "currently only 2D " *
        "HPolygon projection is supported"

    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type, ε)
end

"""
    project(S::LazySet, block::AbstractVector{Int}, ε::Real, [n]::Int=dim(S);
            [kwargs...])

Project a set to a given block and set type with a certified error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `ε`     -- error bound for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set with the given error bound `ε`.
The target set type is chosen automatically.
"""
@inline function project(S::LazySet, block::AbstractVector{Int}, ε::Real,
                         n::Int=dim(S); kwargs...)
    # currently we only support HPolygon
    if length(block) == 2
        set_type = HPolygon
    else
        throw(ArgumentError("ε-close approximation is only supported for 2D " *
                            "blocks"))
    end
    return project(S, block, set_type => ε, n)
end

"""
    rectify(X::LazySet, [concrete_intersection]::Bool=false)

Concrete rectification of a set.

### Input

- `X`                     -- set
- `concrete_intersection` -- (optional, default: `false`) flag to compute
                             concrete intersections for intermediate results

### Output

A set corresponding to the rectification of `X`, which is in general a union of
linear maps of intersections.

### Algorithm

For each dimension in which `X` is both positive and negative, we split `X` into
these two parts. Additionally we project the negative part to zero.
"""
function rectify(X::LazySet, concrete_intersection::Bool=false)
    return to_union_of_projections(Rectification(X), concrete_intersection)
end

"""
    rationalize(::Type{T}, X::LazySet{<:AbstractFloat}, tol::Real)
        where {T<:Integer}

Approximate a set of floating-point numbers as a set whose entries are rationals
of the given integer type.

### Input

- `T`   -- (optional, default: `Int`) integer type to represent the rationals
- `X`   -- set which has floating-point components
- `tol` -- (optional, default: `eps(N)`) tolerance of the result; each rationalized
           component will differ by no more than `tol` with respect to the floating-point value

### Output

A set of the same base type of `X` where each numerical component is of
type `Rational{T}`.
"""
function rationalize(::Type{T}, X::LazySet{<:AbstractFloat}, tol::Real) where {T<:Integer}
    m = length(fieldnames(typeof(X)))
    frat = ntuple(fi -> rationalize(T, getfield(X, fi), tol), m)
    ST = basetype(X)
    return ST(frat...)
end

# no integer type specified
rationalize(X::LazySet{<:AbstractFloat}; kwargs...) =
    rationalize(Int, X; kwargs...)

# `tol` as kwarg
rationalize(::Type{T}, X::LazySet{N}; tol::Real=eps(N)) where {T<:Integer, N<:AbstractFloat} =
    rationalize(T, X, tol)

# vectors of sets
rationalize(::Type{T}, X::AbstractVector{<:LazySet{<:AbstractFloat}}, tol::Real) where {T<:Integer} =
    rationalize.(Ref(T), X, Ref(tol))

"""
    permute(X::LazySet, p::AbstractVector{Int})

Permute the dimensions of a set according to a given permutation vector.

### Input

- `X` -- set
- `p` -- permutation vector

### Output

A new set corresponding to `X` where the dimensions have been permuted according
to `p`.
"""
function permute end

"""
    chebyshev_center_radius(P::LazySet{N};
                            [backend]=default_polyhedra_backend(P),
                            [solver]=default_lp_solver_polyhedra(N; presolve=true)
                           ) where {N}

Compute a [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of a polytopic set.

### Input

- `P`       -- polytopic set
- `backend` -- (optional; default: `default_polyhedra_backend(P)`) the backend
               for polyhedral computations
- `solver`  -- (optional; default:
               `default_lp_solver_polyhedra(N; presolve=true)`) the LP solver
               passed to `Polyhedra`

### Output

The pair `(c, r)` where `c` is a Chebyshev center of `P` and `r` is the radius
of the largest ball with center `c` enclosed by `P`.

### Notes

The Chebyshev center is the center of a largest Euclidean ball enclosed by `P`.
In general, the center of such a ball is not unique, but the radius is.

### Algorithm

We call `Polyhedra.chebyshevcenter`.
"""
function chebyshev_center_radius(P::LazySet{N};
                                 backend=default_polyhedra_backend(P),
                                 solver=default_lp_solver_polyhedra(N; presolve=true)
                                ) where {N}
    require(@__MODULE__, :Polyhedra; fun_name="chebyshev_center")
    if !is_polyhedral(P) && !isboundedtype(typeof(P))
        error("can only compute a Chebyshev center for polytopes")
    end

    Q = polyhedron(P; backend=backend)
    c, r = Polyhedra.chebyshevcenter(Q, solver)
    return c, r
end

function load_polyhedra_lazyset()  # function to be loaded by Requires
return quote
# see the interface file init_Polyhedra.jl for the imports

"""
    polyhedron(P::LazySet; [backend]=default_polyhedra_backend(P))

Compute a set representation from `Polyhedra.jl`.

### Input

- `P`       -- polyhedral set
- `backend` -- (optional, default: call `default_polyhedra_backend(P)`)
                the polyhedral computations backend

### Output

A set representation in the `Polyhedra` library.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).

### Algorithm

This default implementation uses `tosimplehrep`, which computes the constraint
representation of `P`. Set types preferring the vertex representation should
implement their own method.
"""
function polyhedron(P::LazySet; backend=default_polyhedra_backend(P))
    A, b = tosimplehrep(P)
    return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
end

"""
    triangulate(X::LazySet)

Triangulate a three-dimensional polyhedral set.

### Input

- `X` -- three-dimensional polyhedral set

### Output

A tuple `(p, c)` where `p` is a matrix, with each column containing a point, and
`c` is a list of 3-tuples containing the indices of the points in each triangle.
"""
function triangulate(X::LazySet)
    dim(X) == 3 || throw(ArgumentError("the dimension of the set should be " *
        "three, got $(dim(X))"))
    @assert is_polyhedral(X) "triangulation requires a polyhedral set"

    P = polyhedron(X)
    mes = Mesh(P)
    coords = Polyhedra.GeometryBasics.coordinates(mes)
    connec = Polyhedra.GeometryBasics.faces(mes)

    ntriangles = length(connec)
    npoints = length(coords)
    @assert npoints == 3 * ntriangles
    points = Matrix{Float32}(undef, 3, npoints)

    for i in 1:npoints
        points[:, i] .= coords[i].data
    end

    connec_tup = getfield.(connec, :data)

    return points, connec_tup
end

end end  # quote / load_polyhedra_lazyset()
