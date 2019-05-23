for N in [Float64, Rational{Int}, Float32]
    # ==============================
    # Testing box approximation
    # ==============================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], N(0.1))
    h = box_approximation(b)
    hexp = Hyperrectangle(N[1, 1], N[0.1, 0.1])
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

    # empty set
    E = EmptySet{N}()
    @test box_approximation(E) == E

    # rectification
    @test box_approximation(Rectification(EmptySet{N}())) isa EmptySet{N}
    r = Rectification(Ball1(N[0, 0], N(1)))
    @test box_approximation(r) == Hyperrectangle(low=N[0, 0], high=N[1, 1])

    # ===================================================================
    # Testing box_approximation_symmetric (= symmetric interval hull)
    # ===================================================================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], to_N(N, 0.1))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], to_N(N, [1.1, 1.1]))
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 2D unit ball in the 1-norm
    b = Ball1(N[1, -2], to_N(N, 0.2))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], to_N(N, [1.2, 2.2]))
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 3D unit ball in the 1-norm
    b = Ball1(N[1, 2, 0], to_N(N, 0.1))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0, 0], to_N(N, [1.1, 2.1, 0.1]))
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 4D hyperrectangle
    b = Hyperrectangle(to_N(N, [-1.5, -2.5, 2.4, -0.4]),
                       to_N(N, [0.1, 0.2, 0.3, 0.4]))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(zeros(N, 4), to_N(N, [1.6, 2.7, 2.7, 0.8]))
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Testing alias symmetric_interval_hull
    h = symmetric_interval_hull(b)

    # empty set
    E = EmptySet{N}()
    @test box_approximation_symmetric(E) == E
end
