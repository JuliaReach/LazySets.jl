```@contents
Pages = ["AbstractBallp.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Balls in the p-norm (AbstractBallp)](@id def_AbstractBallp)

A ball is a centrally-symmetric set with a characteristic p-norm.

```@docs
AbstractBallp
```

This interface requires to implement the following functions:

```@docs
radius_ball(::AbstractBallp)
ball_norm(::AbstractBallp)
```

This interface defines the following functions:

```@docs
σ(::AbstractVector, ::AbstractBallp)
ρ(::AbstractVector, ::AbstractBallp)
∈(::AbstractVector, ::AbstractBallp)
```

## Implementations

* [Ball1](@ref def_Ball1)
* [Ball2](@ref def_Ball2)
* [BallInf](@ref def_BallInf)
* [Ballp](@ref def_Ballp)
