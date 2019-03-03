import Base: ==, copy

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
       translate

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

```jldoctest
julia> subtypes(LazySet)
17-element Array{Any,1}:
 AbstractCentrallySymmetric
 AbstractPolyhedron
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
"""
abstract type LazySet{N} end


# --- common LazySet functions ---


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
            d = LazySets.Approximations.UnitVector(i, n, o)
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

Return whether two LazySets of the same type are exactly equal by recursively
comparing their fields until a mismatch is found.

### Input

- `X` -- any `LazySet`
- `Y` -- another `LazySet` of the same type as `X`

### Output

- `true` iff `X` is equal to `Y`.

### Notes

The check is purely syntactic and the sets need to have the same base type.
I.e. `X::VPolytope == Y::HPolytope` returns `false` even if `X` and `Y` represent the
same polytope. However `X::HPolytope{Int64} == Y::HPolytope{Float64}` is a valid comparison.

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
    if Compat.isabstracttype(promote_type(typeof(X), typeof(Y)))
        return false
    end

    for f in fieldnames(typeof(X))
        if getfield(X, f) != getfield(Y, f)
            return false
        end
    end

    return true
end

@static if VERSION >= v"0.7-"
    # hook into random API
    import Random.rand
    function rand(rng::AbstractRNG, ::SamplerType{T}) where T<:LazySet
        rand(T, rng=rng)
    end
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
    if X isa Universe
        result = true
    elseif X isa HPolyhedron
        # the polyhedron is universal iff there is no constraint at all
        result = isempty(constraints_list(X))
    elseif X isa HalfSpace || X isa Hyperplane || X isa Line
        result = false
    elseif isbounded(X)
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
