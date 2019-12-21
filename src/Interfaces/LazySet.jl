import Base: ==, ≈, copy, eltype
import Random.rand

export LazySet,
       ρ, support_function,
       σ, support_vector,
       dim,
       norm,
       radius,
       diameter,
       an_element,
       isbounded, isbounded_unit_dimensions,
       neutral,
       absorbing,
       tosimplehrep,
       isuniversal,
       translate,
       affine_map,
       is_interior_point,
       isoperation,
       isoperationtype,
       isequivalent,
       isconvextype

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
- `dim(S::LazySet)::Int` -- the ambient dimension of `S`

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
19-element Array{Any,1}:
 AbstractCentrallySymmetric
 AbstractPolyhedron
 AffineMap
 Bloating
 CacheMinkowskiSum
 CartesianProduct
 CartesianProductArray
 ConvexHull
 ConvexHullArray
 EmptySet
 ExponentialMap
 ExponentialProjectionMap
 Intersection
 IntersectionArray
 LinearMap
 MinkowskiSum
 MinkowskiSumArray
 ResetMap
 Translation
```

If we only consider *concrete* subtypes, then:

```jldoctest; setup = :(using LazySets: subtypes)
julia> concrete_subtypes = subtypes(LazySet, true);

julia> length(concrete_subtypes)
39

julia> println.(concrete_subtypes);
AffineMap
Ball1
Ball2
BallInf
Ballp
Bloating
CacheMinkowskiSum
CartesianProduct
CartesianProductArray
ConvexHull
ConvexHullArray
Ellipsoid
EmptySet
ExponentialMap
ExponentialProjectionMap
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
Line
LineSegment
LinearMap
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
    ρ(d::AbstractVector{N}, S::LazySet{N})::N where {N<:Real}

Evaluate the support function of a set in a given direction.

### Input

- `d` -- direction
- `S` -- convex set

### Output

The support function of the set `S` for the direction `d`.

### Notes

The numeric type of the direction and the set must be identical.
"""
function ρ(d::AbstractVector{N}, S::LazySet{N})::N where {N<:Real}
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
    isbounded(S::LazySet)::Bool

Determine whether a set is bounded.

### Input

- `S` -- set

### Output

`true` iff the set is bounded.

### Algorithm

We check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(S::LazySet)::Bool
    return isbounded_unit_dimensions(S)
end

"""
    isbounded_unit_dimensions(S::LazySet{N})::Bool where {N<:Real}

Determine whether a set is bounded in each unit dimension.

### Input

- `S` -- set

### Output

`true` iff the set is bounded in each unit dimension.

### Algorithm

This function performs ``2n`` support function checks, where ``n`` is the
ambient dimension of `S`.
"""
function isbounded_unit_dimensions(S::LazySet{N})::Bool where {N<:Real}
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
        return radius(Approximations.ballinf_approximation(S)::BallInf, p)
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
    affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector)

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
function affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector)
    return translate(linear_map(M, X), v)
end

"""
    an_element(S::LazySet{N}) where {N<:Real}

Return some element of a convex set.

### Input

- `S` -- convex set

### Output

An element of a convex set.
"""
function an_element(S::LazySet{N}) where {N<:Real}
    return σ(sparsevec([1], [one(N)], dim(S)), S)
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
    isuniversal(X::LazySet{N}, [witness]::Bool=false
               )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function isuniversal(X::LazySet{N}, witness::Bool=false
                    )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
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

### Notes

Plotting of unbounded sets is not implemented yet (see
[#576](https://github.com/JuliaReach/LazySets.jl/issues/576)).

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
