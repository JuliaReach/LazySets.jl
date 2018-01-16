# A Reachability Algorithm Using Zonotopes

```@contents
Pages = ["reach_zonotopes.md"]
Depth = 3
```

## Introduction

In this section we present an algorithm implemented using `LazySets` that computes
the reach sets of an affine ordinary differential equation (ODE). This
algorithm is from A. Girard's *"Reachability of uncertain linear systems using zonotopes*,
HSCC. Vol. 5. 2005. We have chosen this algorithm for the purpose of illustration
of a complete application of `LazySets`.

Let us introduce some notation. Consider the continuous initial set-valued problem
(IVP)

```math
    x'(t) = A x(t) + u(t)
```
in the time interval ``t ∈ [0, T]``, where:

-  ``A`` is a real matrix of order ``n``,
- ``u(t)`` is a non-deterministic input such that ``\Vert u(t) \Vert_∞ ≦ μ`` for all ``t``,
- ``x(0) ∈ \mathcal{X}_0``, where ``\mathcal{X}_0`` is a convex set.

Given a step size ``δ``, `Algorithm1` returns a sequence of sets that overapproximates
the states reachable by any trajectory of this IVP.

## Algorithm

```@example example_reach_zonotopes
using LazySets, Plots

function Algorithm1(A, X0, δ, μ, T)

    # bloating factors
    Anorm = norm(A, Inf)
    α = (expm(δ*Anorm) - 1 - δ*Anorm)/norm(X0, Inf)
    β = (expm(δ*Anorm) - 1)*μ/Anorm

    # discretized system
    n = size(A, 1)
    ϕ = expm(δ*A)
    N = floor(Int, T/δ)

    # preallocate arrays
    Q = Vector{LazySet}(N)
    R = Vector{LazySet}(N)

    # initial reach set in the time interval [0, δ]
    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    c = X0.center
    Q1_generators = hcat(ϕp * X0.generators, ϕm * c, ϕm * X0.generators)
    Q[1] = Zonotope(ϕp * c, Q1_generators) ⊕ BallInf(zeros(n), α + β)
    R[1] = Q[1]

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = BallInf(zeros(n), β)
    for i in 2:N
        Q[i] = ϕ * Q[i-1] ⊕ ballβ
        R[i] = Q[i]
    end
    return R
end
```

## Projection

```@example example_reach_zonotopes
function project(R, vars, n)
    # projection matrix
    M = sparse(1:2, vars, [1., 1.], 2, n)
    return [M * Ri for Ri in R]
end
```

## Example 1

```@example example_reach_zonotopes
A = [-1 -4; 4 -1]
X0 = Zonotope([1.0, 0.0], 0.1*eye(2))
μ = 0.05
δ = 0.02
T = 2.

R = Algorithm1(A, X0, δ, μ, 2.*δ); # warm-up

R = Algorithm1(A, X0, δ, μ, T)

plot(R, 1e-2, fillalpha=0.1)
```


## Example 2

```@example example_reach_zonotopes
A = Matrix{Float64}([-1 -4 0 0 0;
                      4 -1 0 0 0;
                      0 0 -3 1 0;
                      0 0 -1 -3 0;
                      0 0 0 0 -2])
X0 = Zonotope([1.0, 0.0, 0.0, 0.0, 0.0], 0.1*eye(5))
μ = 0.01
δ = 0.005
T = 1.

R = Algorithm1(A, X0, δ, μ, 2*δ); # warm-up

R = Algorithm1(A, X0, δ, μ, T)
R = project(R, [1, 3], 5)

plot(R, 1e-2, fillalpha=0.1, xlabel="x1", ylabel="x3")
```
