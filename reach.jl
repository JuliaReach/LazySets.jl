using LazySets

function reach_hybrid(As, Gs, init, δ, μ, T, max_order)
    queue = [(init[1], init[2], 0)]

    res = Zonotope[]
    while !isempty(queue)
        init, loc, t = pop!(queue)
        R = reach_continuous_ordred(As[loc], init, δ, μ, T-t, max_order)
        append!(res, R)
        for (guard, tgt_loc) in Gs[loc]
            for i in 1:length(R)
                S = R[i]
                if !is_intersection_empty(S, guard)
                    push!(queue, (S, tgt_loc, δ * (i-1)))
                end
            end
        end
    end
    return res
end



function reach_continuous(A, X0, δ, μ, T)

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

function reach_continuous_ordred(A, X0, δ, μ, T, max_order)

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
    Q[1] = minkowski_sum(Zonotope(ϕp * c, Q1_generators), Zonotope(zeros(n), (α + β)*eye(n)))
    R[1] = Q[1]
    init_order = order(Q[1])

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = Zonotope(zeros(n), β*eye(n))
    for i in 2:N
        Q[i] = minkowski_sum(linear_map(ϕ, Q[i-1]), ballβ)
        if order(Q[i]) > max_order
            Q[i] = reduce_order(Q[i], init_order)
        end
        R[i] = Q[i]
    end
    return R
end

function example()
    # dynamics
    A1 = A2 = [-1.0 0.0; 1.0 0.0]
    A3 = A4 = [-1.0 0.0; 1.0 -1.0]
    As = [A1, A2, A3, A4]

    # transitions
    t1 = [(Hyperplane([1., 0.], -1.), 2), (Hyperplane([0., 1.], 1.), 3)]
    t2 = [(Hyperplane([0., 1.], 1.), 3)]
    t3 = [(Hyperplane([0., 1.], 0.), 1), (Hyperplane([1., 0.], -1.), 4)]
    t4 = [(Hyperplane([0., 1.], 0.),2), (Hyperplane([1., 0.], 1.), 3)]
    Gs = [t1, t2, t3, t4]

    # initial condition
    X0 = Zonotope([2., 1.], [[0.5, 0.]])
    init_loc = 3
    init = (X0, init_loc)

    # input uncertainty
    μ = 0.05

    # discretization step
    δ = 0.02

    # time bound
    T = 2.

    # run analysis
    reach_hybrid(As, Gs, init, δ, μ, T, 40)
end
