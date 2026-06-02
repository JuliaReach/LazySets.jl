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

# helper function to compute the bloating ball
function _bloating_ball(B::Bloating)
    return _bloating_ball(B.ε, B.p, dim(B))
end

function _bloating_ball(ε::N, p::N, n::Int) where {N}
    @assert ε >= zero(N) "cannot compute the ball for a negative bloating"
    return Ballp(p, zeros(N, n), ε)
end

include("an_element.jl")
include("center.jl")
include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
