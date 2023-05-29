# A Hybrid Reachability Algorithm Using Zonotopes

```@contents
Pages = ["reach_zonotopes_hybrid.md"]
Depth = 3
```

## Introduction

In this section we present an algorithm implemented using `LazySets` that
computes the reach sets of a hybrid system of linear ordinary differential
equations (ODE).
This algorithm is an extension of the one presented in
[A Reachability Algorithm Using Zonotopes](@ref).

We consider a simple case here where modes do not have invariants and
transitions do not have updates.
In set-based analysis like ours, it may make sense to take a transition as soon
as one state in the current set of states can take it.
Note that this is not equivalent to *must* semantics of hybrid automata (also
called *urgent transitions*), which is defined on single trajectories.
We also offer the usual *may* transitions interpretation.


## Hybrid algorithm

The hybrid algorithm maintains a queue of triples ``(m, X, t)`` where ``m`` is a
mode, ``X`` is a set of states, and ``t`` is a time point.
For each element in the queue the algorithm calls the
[Continuous algorithm](@ref) to compute the reachable states in the current mode
``m``, starting in the current states ``X`` at time ``t``.
The result is a flowpipe, i.e., a sequence of sets of states.
For each of those sets we check intersection with the guards of ``m``'s outgoing
transitions.
Depending on the transition semantics, we add the discrete successors to the
queue and continue with the next iteration until the queue is empty.

```@example example_reach_zonotopes_hybrid
using Plots, LazySets, LinearAlgebra

function reach_hybrid(As, Ts, init, δ, μ, T, max_order, instant_transitions)
    # initialize queue with initial mode and states at time t=0
    queue = Vector{Tuple{LazySet, Integer, Float64}}(undef, 1)
    queue[1] = (init[1], init[2], 0.0)

    res = Tuple{LazySet, Int}[]
    while !isempty(queue)
        init, loc, t = pop!(queue)
        println("currently in location $loc at time $t")
        R = reach_continuous(As[loc], init, δ, μ, T-t, max_order)
        found_transition = false
        for i in 1:length(R)-1
            S = R[i]
            push!(res, (S, loc))
            for (guard, tgt_loc) in Ts[loc]
                if !isdisjoint(S, guard)
                    new_t = t + δ * i
                    push!(queue, (S, tgt_loc, new_t))
                    found_transition = true
                    println("transition $loc -> $tgt_loc at time $new_t")
                end
            end
            if instant_transitions && found_transition
                break
            end
        end
        if !instant_transitions || !found_transition && length(R) > 0
            push!(res, (R[end], loc))
        end
    end
    return res
end
nothing # hide
```

### Continuous algorithm

This is basically the same implementation as outlined in the section
[A Reachability Algorithm Using Zonotopes](@ref), only that this time we use
concrete operations on zonotopes.

```@example example_reach_zonotopes_hybrid
function reach_continuous(A, X0, δ, μ, T, max_order)
    # bloating factors
    Anorm = norm(A, Inf)
    α = (exp(δ*Anorm) - 1 - δ*Anorm)/norm(X0, Inf)
    β = (exp(δ*Anorm) - 1)*μ/Anorm

    # discretized system
    n = size(A, 1)
    ϕ = exp(δ*A)
    N = floor(Int, T/δ)

    # preallocate array
    R = Vector{Zonotope}(undef, N)
    if N == 0
        return R
    end

    # initial reach set in the time interval [0, δ]
    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    c = X0.center
    gens = hcat(ϕp * X0.generators, ϕm * c, ϕm * X0.generators)
    R[1] = minkowski_sum(Zonotope(ϕp * c, gens),
                         Zonotope(zeros(n), Matrix((α + β)*I, n, n)))
    if order(R[1]) > max_order
        R[1] = reduce_order(R[1], max_order)
    end

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = Zonotope(zeros(n), Matrix(β*I, n, n))
    for i in 2:N
        R[i] = minkowski_sum(linear_map(ϕ, R[i-1]), ballβ)
        if order(R[i]) > max_order
            R[i] = reduce_order(R[i], max_order)
        end
    end
    return R
end
nothing # hide
```

### Plotting results

For illustration purposes it is helpful to plot the flowpipes in different
colors, depending on the current mode.
The following function does that for 2-mode models.

```@example example_reach_zonotopes_hybrid
function plot_res(res)
    p = plot()
    for i in 1:length(res)
        if res[i][2] == 1
            c = "blue"
        elseif res[i][2] == 2
            c = "red"
        end
        plot!(p, reduce_order(res[i][1], 2), color=c, alpha=0.1)
    end
    return p
end
nothing # hide
```

## Example

We consider an extension of the example presented in
*Reachability of uncertain linear systems using zonotopes*, A. Girard, HSCC.
Vol. 5. 2005 to a hybrid system with two modes ``\ell_i``, ``i = 1, 2``, with
initial states
``[0.9, 1.1] \times [-0.1, 0.1]`` and uncertain inputs from a set ``u`` with
``\mu = \Vert u \Vert_\infty = 0.001``.

The dynamics matrices ``A_i`` are defined as follows:

```math
	A_1 = \begin{pmatrix} -1 & -4 \\ 4 & -1 \end{pmatrix}, \qquad A_2 = \begin{pmatrix} 1 & 4 \\ -4 & -1 \end{pmatrix}.
```
We add a transition ``t_i`` from mode ``\ell_i`` to ``\ell_{3-i}`` with a
hyperplane guard ``g_i``:

```math
	g_1 \triangleq x_1 = -0.5 \qquad g_2 \triangleq x_2 = -0.3
```

`LazySets` offers an order reduction function for zonotopes, which we used here
with an upper bound of 10 generators. We plot the reachable states for the time
interval ``[0, 4]`` and time step ``δ = 0.01``.

```@example example_reach_zonotopes_hybrid
    # dynamics
    A1 = [-1 -4; 4 -1]
    A2 = [1 4; -4 -1]
    As = [A1, A2]

    # transitions
    t1 = [(Hyperplane([1., 0.], -0.5), 2)]
    t2 = [(Hyperplane([0., 1.], -0.3), 1)]
    Ts = [t1, t2]

    # initial condition
    X0 = Zonotope([1.0, 0.0], Matrix(0.1*I, 2, 2))
    init_loc = 1
    init = (X0, init_loc)

    # input uncertainty
    μ = 0.001

    # discretization step
    δ = 0.01

    # time bound
    T = 4.

    # maximum order of zonotopes
    max_order = 10

    # take transitions only the first time they are enabled?
    instant_transitions = true

    # run analysis
    res = reach_hybrid(As, Ts, init, δ, μ, T, max_order, instant_transitions)

    # plot result
    plot_res(res)
```
