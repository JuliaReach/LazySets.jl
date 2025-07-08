for N in [Float64, Rational{Int}, Float32]
    # ==============================
    # Testing box approximation
    # ==============================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], N(1 // 10))
    h = box_approximation(b)
    hexp = Hyperrectangle(N[1, 1], N[1 // 10, 1 // 10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 2D unit ball in the 1-norm
    b = Ball1(N[1, -2], N(0.2))
    h = box_approximation(b)
    hexp = Hyperrectangle(N[1, -2], N[0.2, 0.2])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 3D unit ball in the 1-norm
    b = Ball1(N[1, 2, 0], N(1))
    h = box_approximation(b)
    hexp = Hyperrectangle(N[1, 2, 0], N[1, 1, 1])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a Hyperrectangle
    hexp = Hyperrectangle(N[0, 0], N[1, 1])
    h = box_approximation(hexp)
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a VPolytope
    P = VPolytope([N[0, 0], N[2, 1], N[1, 2]])
    H = box_approximation(P)
    @test H == Hyperrectangle(N[1, 1], N[1, 1])

    # empty set
    E = EmptySet{N}(2)
    @test box_approximation(E) == E

    # alias constructors
    @test interval_hull(b) == □(b) == box_approximation(b)

    # box_approximation of lazy intersection
    # both arguments are bounded
    X = Ball1(zeros(N, 2), N(1))
    Y = Ball1(N[1 // 2, 0], N(1))
    Z = box_approximation(X ∩ Y)
    @test isequivalent(Z, Hyperrectangle(N[1 // 4, 0], N[3 // 4, 3 // 4]))
    # one argument is unbounded
    X = Hyperrectangle(N[0], N[1]) × Universe{N}(1)
    Y = Hyperrectangle(N[0, 0], N[1, 1])
    Z = box_approximation(X ∩ Y)
    @test isequivalent(Z, Y)

    # ===================================================================
    # Testing box_approximation_symmetric (= symmetric interval hull)
    # ===================================================================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], N(1 // 10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], N[11 // 10, 11 // 10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 2D unit ball in the 1-norm
    b = Ball1(N[1, -2], N(2 // 10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], N[12 // 10, 22 // 10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 3D unit ball in the 1-norm
    b = Ball1(N[1, 2, 0], N(1 // 10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0, 0], N[11 // 10, 21 // 10, 1 // 10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 4D hyperrectangle
    b = Hyperrectangle(N[-15 // 10, -25 // 10, 24 // 10, -4 // 10],
                       N[1 // 10, 2 // 10, 3 // 10, 4 // 10])
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(zeros(N, 4), N[16 // 10, 27 // 10, 27 // 10, 8 // 10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Testing alias symmetric_interval_hull
    h == symmetric_interval_hull(b)

    # empty set
    E = EmptySet{N}(2)
    @test box_approximation_symmetric(E) == E
end

for N in [Float64]
    # slightly contradicting bounds are interpreted as a flat set
    X = HalfSpace(N[-1], N(0)) ∩ HalfSpace(N[1], N(-1e-15))
    @test box_approximation(X) == Hyperrectangle(N[-5.0e-16], N[0])

    # box approximation of Taylor model
    # (currently gives different result for non-Float64:
    #  https://github.com/JuliaIntervals/TaylorModels.jl/issues/158)
    I = IA.interval(N(0), N(0))  # interval remainder
    # TaylorModel1
    t = TaylorModels.Taylor1(3)
    q₁ = 1 + 2 * t + 2 * t^2
    D = IA.interval(N(-1), N(1))
    local x0 = IA.mid(D)
    local vTM = [TaylorModels.TaylorModel1(q₁, I, x0, D)]
    @test box_approximation(vTM) == Hyperrectangle(N[2], N[3])
    # TaylorModelN
    local x₁, x₂, x₃ = TaylorModels.set_variables(N, ["x₁", "x₂", "x₃"]; order=5)
    local p₁ = 1 + x₁ - x₂
    local p₂ = x₃ - x₁
    Dx₁ = IA.interval(N(-1), N(1))
    Dx₂ = IA.interval(N(-1), N(1))
    Dx₃ = IA.interval(N(-1), N(1))
    D = Dx₁ × Dx₂ × Dx₃
    local x0 = IntervalBox(IA.mid.(D)...)
    local vTM = [TaylorModels.TaylorModelN(pi, I, x0, D) for pi in [p₁, p₂]]
    @test box_approximation(vTM) == Hyperrectangle(N[1, 0], N[2, 2])
end
