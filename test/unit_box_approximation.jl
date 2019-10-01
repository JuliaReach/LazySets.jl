for N in [Float64, Rational{Int}, Float32]
    # ==============================
    # Testing box approximation
    # ==============================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], N(1//10))
    h = box_approximation(b)
    hexp = Hyperrectangle(N[1, 1], N[1//10, 1//10])
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

    # ===================================================================
    # Testing box_approximation_symmetric (= symmetric interval hull)
    # ===================================================================
    # Box approximation of a 2D square
    b = BallInf(N[1, 1], N(1//10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], N[11//10, 11//10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 2D unit ball in the 1-norm
    b = Ball1(N[1, -2], N(2//10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0], N[12//10, 22//10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 3D unit ball in the 1-norm
    b = Ball1(N[1, 2, 0], N(1//10))
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(N[0, 0, 0], N[11//10, 21//10, 1//10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Box approximation of a 4D hyperrectangle
    b = Hyperrectangle(N[-15//10, -25//10, 24//10, -4//10],
                       N[1//10, 2//10, 3//10, 4//10])
    h = box_approximation_symmetric(b)
    hexp = Hyperrectangle(zeros(N, 4), N[16//10, 27//10, 27//10, 8//10])
    @test h.center ≈ hexp.center
    @test h.radius ≈ hexp.radius

    # Testing alias symmetric_interval_hull
    h = symmetric_interval_hull(b)

    # empty set
    E = EmptySet{N}()
    @test box_approximation_symmetric(E) == E
end

for N in [Float64, Float32]
    # empty intersection due to line search
    X = HalfSpace(N[-1], N(0)) ∩ HalfSpace(N[1], N(-1e-15))
    @test box_approximation(X) isa EmptySet{N}
end
