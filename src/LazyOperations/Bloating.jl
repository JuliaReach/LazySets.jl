export Bloating

"""
    Bloating{N, S<:LazySet{N}} <:LazySet{N}

Type that represents a uniform expansion of a set in a given norm (also
known as *bloating*).

### Fields

- `X` -- set
- `ε` -- (positive) bloating factor
- `p` -- ``p``-norm (should be ``≥ 1``; default: ``2``)

### Notes

`Bloating(X, ε, p)` is equivalent to the Minkowski sum of `X` and a ball in the
`p`-norm of radius `ε` centered in the origin `O` (i.e., `X ⊕ Ballp(p, O, ε)`).
"""
struct Bloating{N, S<:LazySet{N}} <:LazySet{N}
    X::S
    ε::N
    p::N

    function Bloating(X::S, ε::N, p::N=N(2)) where {N, S<:LazySet{N}}
        @assert ε > zero(N) "bloating requires a distance > 0, but $ε was given"
        @assert p >= one(N) "bloating requires a norm >= 1, but $p was given"

        return new{N, S}(X, ε, p)
    end
end

isoperationtype(::Type{<:Bloating}) = true
isconvextype(::Type{Bloating{N, S}}) where {N, S} = isconvextype(S)

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
function _bloating_ball(B::Bloating{N}) where {N}
    return Ballp(B.p, zeros(N, dim(B)), B.ε)
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
function σ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support vector in the zero direction is undefined"

    return σ(d, B.X) + σ(d, _bloating_ball(B))
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
function ρ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support function in the zero direction is undefined"

    return ρ(d, B.X) + ρ(d, _bloating_ball(B))
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

"""
    isempty(B::Bloating)

Determine whether a bloated set is empty.

### Input

- `B` -- bloated set

### Output

`true` iff the wrapped set is empty.
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

The implementation returns the result of `an_element` for the wrapped set.
"""
function an_element(B::Bloating)
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

The constraints list is only available for bloating in the `p`-norm for
``p = 1`` or ``p = ∞`` and if `constraints_list` is available for the unbloated
set.

### Algorithm

We call `constraints_list` on the lazy Minkowski sum.
"""
function constraints_list(B::Bloating)
    @assert (B.p == 1 || B.p == Inf) "the constraints list is only available " *
        "for bloating in the 1-norm or in the infinity norm"

    return constraints_list(MinkowskiSum(B.X, _bloating_ball(B)))
end
