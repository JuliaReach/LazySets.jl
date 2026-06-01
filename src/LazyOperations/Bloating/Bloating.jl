export Bloating

"""
    Bloating{N, S<:LazySet{N}} <: LazySet{N}

Type that represents a uniform expansion of a set in a given norm (also
known as *bloating*).

### Fields

- `X` -- set
- `ε` -- (usually positive) bloating factor
- `p` -- ``p``-norm (should be ``≥ 1``; default: ``2``)

### Notes

The `Bloating` operation preserves convexity: if `X` is convex, then any
bloating of `X` is convex as well.

If `ε` is positive, then `Bloating(X, ε, p)` is equivalent to the Minkowski sum
of `X` and a ball in the `p`-norm of radius `ε` centered in the origin `O`
(i.e., `X ⊕ Ballp(p, O, ε)`).

Some operations require, or silently assume, that `ε` is positive. Check the
documentation for further information.
"""
struct Bloating{N,S<:LazySet{N}} <: LazySet{N}
    X::S
    ε::N
    p::N

    function Bloating(X::S, ε::N, p::N=N(2)) where {N,S<:LazySet{N}}
        @assert p >= one(N) "bloating requires a norm >= 1, but $p was given"

        return new{N,S}(X, ε, p)
    end
end

isoperationtype(::Type{<:Bloating}) = true
isconvextype(::Type{Bloating{N,S}}) where {N,S} = isconvextype(S)

function concretize(B::Bloating{N}) where {N}
    return minkowski_sum(concretize(B.X), Ballp(B.p, zeros(N, dim(B)), B.ε))
end

"""
    dim(B::Bloating)

Return the dimension of a bloated set.

### Input

- `B` -- bloated set

### Output

The ambient dimension of the bloated set.
"""
function dim(B::Bloating)
    return dim(B.X)
end

# helper function to compute the bloating ball
function _bloating_ball(B::Bloating)
    return _bloating_ball(B.ε, B.p, dim(B))
end

function _bloating_ball(ε::N, p::N, n::Int) where {N}
    @assert ε >= zero(N) "cannot compute the ball for a negative bloating"
    return Ballp(p, zeros(N, n), ε)
end

"""
    σ(d::AbstractVector, B::Bloating)

Return the support vector of a bloated set in a given direction.

### Input

- `d` -- direction
- `B` -- bloated set

### Output

The support vector of the bloated set in the given direction.
"""
@validate function σ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support vector in the zero direction is undefined"
    @assert B.ε >= 0 || B.p > 1 "the support vector for negative bloating " *
                                "in the 1-norm is not implemented"

    return σ(d, B.X) +
           sign_cadlag(B.ε) * σ(d, _bloating_ball(abs(B.ε), B.p, dim(B)))
end

"""
    ρ(d::AbstractVector, B::Bloating)

Return the support function of a bloated set in a given direction.

### Input

- `d` -- direction
- `B` -- bloated set

### Output

The support function of the bloated set in the given direction.
"""
@validate function ρ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support function in the zero direction is undefined"

    return ρ(d, B.X) +
           sign_cadlag(B.ε) * ρ(d, _bloating_ball(abs(B.ε), B.p, dim(B)))
end

"""
    isbounded(B::Bloating)

Determine whether a bloated set is bounded.

### Input

- `B` -- bloated set

### Output

`true` iff the wrapped set is bounded.
"""
function isbounded(B::Bloating)
    return isbounded(B.X)
end

function isboundedtype(::Type{<:Bloating{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
    isempty(B::Bloating)

Determine whether a bloated set is empty.

### Input

- `B` -- bloated set

### Output

`true` iff the wrapped set is empty.

### Notes

This implementation disregards negative bloating, which could potentially turn a
non-empty set into an empty set.
"""
function isempty(B::Bloating)
    return isempty(B.X)
end

"""
    an_element(B::Bloating)

Return some element of a bloated set.

### Input

- `B` -- bloated set

### Output

An element in the bloated set.

### Algorithm

This implementation disregards negative bloating and returns the result of
`an_element` for the wrapped set.
"""
function an_element(B::Bloating)
    if B.ε < 0
        throw(ArgumentError("negative bloating is not supported"))
    end
    return an_element(B.X)
end

"""
    constraints_list(B::Bloating)

Return the list of constraints of a bloated set.

### Input

- `B` -- bloated set

### Output

The list of constraints of the bloated set.

### Notes

The constraints list is only available for non-negative bloating in the `p`-norm
for ``p = 1`` or ``p = ∞`` and if `constraints_list` is available for the
unbloated set.

### Algorithm

We call `constraints_list` on the lazy Minkowski sum with the bloating ball.
"""
function constraints_list(B::Bloating)
    @assert ispolyhedral(B) "the constraints list is only available for " *
                            "polyhedral bloating (which requires a polyhedral base set and the " *
                            "1-norm or the infinity norm)"
    if B.ε < 0
        throw(ArgumentError("computing the constraints list of a negatively " *
                            "bloated set is not supported"))
    end

    return constraints_list(MinkowskiSum(B.X, _bloating_ball(B)))
end

"""
    center(B::Bloating)

Return the center of a bloated set.

### Input

- `B` -- bloated set

### Output

The center of the wrapped set.

### Notes

This implementation disregards negative bloating, which could potentially remove
the center from the set.
"""
function center(B::Bloating)
    return center(B.X)
end

"""
    ispolyhedral(B::Bloating)

Check whether a bloated set is polyhedral.

### Input

- `B` -- bloated set

### Output

`true` if the set is polyhedral.

### Algorithm

We check the sufficient condition that the base set is polyhedral and that the
norm for bloating is either 1-norm or the infinity norm.
"""
function ispolyhedral(B::Bloating)
    return (B.p == 1 || B.p == Inf) && ispolyhedral(B.X)
end

@validate function translate(B::Bloating, x::AbstractVector)
    return Bloating(translate(B.X, x), B.ε, B.p)
end
