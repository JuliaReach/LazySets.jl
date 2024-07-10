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
ball_norm(::AbstractBallp)
radius_ball(::AbstractBallp)
```

This interface defines the following functions:

```@docs
∈(::AbstractVector, ::AbstractBallp)
ρ(::AbstractVector, ::AbstractBallp)
σ(::AbstractVector, ::AbstractBallp)
```

## Implementations

* [Ball1](@ref def_Ball1)
* [Ball2](@ref def_Ball2)
* [BallInf](@ref def_BallInf)
* [Ballp](@ref def_Ballp)
