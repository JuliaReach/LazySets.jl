using LazySets

function Phi1(A, δ)
    n = size(A, 1)
    P = expm(full([A * δ sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)]))
    Phi1Adelta = P[1:n, (n+1):2*n]
end

function reach_hybrid(As, bs, Gs, init, δ, μ, T, max_order, must_semantics)
    queue = [(init[1], init[2], 0.)]

    res = Zonotope[]
    while !isempty(queue)
        init, loc, t = pop!(queue)
        R = reach_continuous(As[loc], bs[loc], init, δ, μ, T-t, max_order)
        found_transition = false
        for i in 1:length(R)-1
            S = R[i]
            push!(res, S)
            for (guard, tgt_loc) in Gs[loc]
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
        if !must_semantics || !found_transition
            push!(res, R[end])
        end
    end
    return res
end

function reach_continuous(A, b, X0, δ, μ, T)

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

function reach_continuous(A, b, X0, δ, μ, T, max_order)

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
    P1 = Phi1(A, δ)
    Q[1] = minkowski_sum(Zonotope(ϕp * c, Q1_generators), Zonotope(P1 * b * δ, (α + β)*eye(n)))
    if order(Q[1]) > max_order
        Q[1] = reduce_order(Q[1], max_order)
    end
    R[1] = Q[1]
#     init_order = order(Q[1])

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = Zonotope(zeros(n), β*eye(n))
    for i in 2:N
        Q[i] = minkowski_sum(minkowski_sum(linear_map(ϕ, Q[i-1]), ballβ), Zonotope(P1 * b * δ, zeros(n, 0)))
        if order(Q[i]) > max_order
            Q[i] = reduce_order(Q[i], max_order)
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
    b1 = [-2.0, 0.0]
    b2 = [3.0, 0.0]
    b3 = [-2.0, -5.0]
    b4 = [3.0, -5.0]
    bs = [b1, b2, b3, b4]

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
    δ = 0.08

    # time bound
    T = 4.

    # maximum order of zonotopes
    max_order = 4

    # use must semantics?
    must_semantics = true

    # run analysis
    reach_hybrid(As, bs, Gs, init, δ, μ, T, max_order, must_semantics)
end
