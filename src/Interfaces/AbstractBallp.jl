export AbstractBallp

"""
    AbstractBallp{N} <: AbstractCentrallySymmetric{N}

Abstract type for p-norm balls.

### Notes

See [`Ballp`](@ref) for a standard implementation of this interface.

Every concrete `AbstractBallp` must define the following methods:

- `center(::AbstractBallp)` -- return the center
- `radius_ball(::AbstractBallp)` -- return the ball radius
- `ball_norm(::AbstractBallp)` -- return the characteristic norm

The subtypes of `AbstractBallp`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractBallp)
2-element Vector{Any}:
 Ball2
 Ballp
```

There are two further set types implementing the `AbstractBallp` interface, but
they also implement other interfaces and hence cannot be subtypes: `Ball1` and
`BallInf`.
"""
abstract type AbstractBallp{N} <: AbstractCentrallySymmetric{N} end
