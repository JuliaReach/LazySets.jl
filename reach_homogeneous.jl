using LazySets

function reach_hybrid(As, Ts, init, δ, μ, T, max_order, must_semantics)
    queue = [(init[1], init[2], 0.)]

    res = Tuple{Zonotope, Int}[]
    while !isempty(queue)
        init, loc, t = pop!(queue)
        println("currently in location $loc at time $t")
        R = reach_continuous(As[loc], init, δ, μ, T-t, max_order)
        found_transition = false
        for i in 1:length(R)-1
            S = R[i]
            push!(res, (S, loc))
            for (guard, tgt_loc) in Ts[loc]
                if !is_intersection_empty(S, guard)
                    new_t = t + δ * i
                    push!(queue, (S, tgt_loc, new_t))
                    found_transition = true
                    println("transition $loc -> $tgt_loc at time $new_t")
                end
            end
            if must_semantics && found_transition
                break
            end
        end
        if (!must_semantics || !found_transition) && length(R) > 0
            push!(res, (R[end], loc))
        end
    end
    return res
end

function reach_continuous(A, X0, δ, μ, T, max_order)
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
    if N == 0
        return R
    end

    # initial reach set in the time interval [0, δ]
    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    c = X0.center
    Q1_generators = hcat(ϕp * X0.generators, ϕm * c, ϕm * X0.generators)
    Q[1] = minkowski_sum(Zonotope(ϕp * c, Q1_generators), Zonotope(zeros(n), (α + β)*eye(n)))
    if order(Q[1]) > max_order
        Q[1] = reduce_order(Q[1], max_order)
    end
    R[1] = Q[1]

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = Zonotope(zeros(n), β*eye(n))
    for i in 2:N
        Q[i] = minkowski_sum(linear_map(ϕ, Q[i-1]), ballβ)
        if order(Q[i]) > max_order
            Q[i] = reduce_order(Q[i], max_order)
        end
        R[i] = Q[i]
    end
    return R
end

function example()
    # dynamics
    A1 = [-1 -4; 4 -1]
    A2 = [1 4; -4 -1]
    As = [A1, A2]

    # transitions
    t1 = [(Hyperplane([1., 0.], -0.5), 2)]
    t2 = [(Hyperplane([0., 1.], -0.3), 1)]
    Ts = [t1, t2]

    # initial condition
    X0 = Zonotope([1.0, 0.0], 0.1*eye(2))
    init_loc = 1
    init = (X0, init_loc)

    # input uncertainty
    μ = 0.001

    # discretization step
    δ = 0.001

    # time bound
    T = 4.

    # maximum order of zonotopes
    max_order = 10

    # use must semantics?
    must_semantics = true

    # run analysis
    reach_hybrid(As, Ts, init, δ, μ, T, max_order, must_semantics)
end

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
