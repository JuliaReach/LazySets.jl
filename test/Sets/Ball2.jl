for N in [Float64, Float32]
    # random ball
    rand(Ball2)

    # 1D Ball2
    b = Ball2(N[0], N(1))
    # Test Dimension
    @test dim(b) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, b) == N[1]
    d = N[-1]
    @test σ(d, b) == N[-1]

    # 2D Ball2
    b = Ball2(N[0, 0], N(1))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, b) == N[1, 0]
    d = N[-1, 0]
    @test σ(d, b) == N[-1, 0]
    d = N[0, 1]
    @test σ(d, b) == N[0, 1]
    d = N[0, -1]
    @test σ(d, b) == N[0, -1]
    d = N[0, 0]
    @test σ(d, b) == N[0, 0]

    # 2D Ball2 not 0-centered
    b = Ball2(N[1, 2], N(1))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, b) == N[2, 2]
    d = N[-1, 0]
    @test σ(d, b) == N[0, 2]
    d = N[0, 1]
    @test σ(d, b) == N[1, 3]
    d = N[0, -1]
    @test σ(d, b) == N[1, 1]
    d = N[0, 0]
    @test σ(d, b) ∈ b

    # 2D Ball2 radius =/= 1
    b = Ball2(N[0, 0], N(2))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, b) == N[2, 0]
    d = N[-1, 0]
    @test σ(d, b) == N[-2, 0]
    d = N[0, 1]
    @test σ(d, b) == N[0, 2]
    d = N[0, -1]
    @test σ(d, b) == N[0, -2]

    # unicode constructor
    @test ○(center(b), radius(b)) == b

    # center
    @test center(b) == N[0, 0]

    # radius_ball
    @test LazySets.radius_ball(b) == N(2)

    # ball_norm
    @test LazySets.ball_norm(b) == N(2)

    # boundedness
    @test isbounded(b) && isboundedtype(typeof(b))

    # ispolyhedral
    @test !ispolyhedral(b)

    # isempty
    @test !isempty(b)

    # isuniversal
    answer, w = isuniversal(b, true)
    @test !isuniversal(b) && !answer && w ∉ b

    # an_element function
    b = Ball2(N[1, 2], N(2))
    @test an_element(b) ∈ b

    # translation
    @test translate(b, N[1, 2]) == Ball2(N[2, 4], N(2))
    bb = Ball2(N[0, 0], N(1))
    @test translate!(bb, N[1, 1]) == Ball2(N[1, 1], N(1)) == bb

    # inclusion of a Ball2 in a Ball2
    b1 = Ball2(N[1, 2], N(2))
    b2 = Ball2(N[1, 2], N(0))
    b3 = Ball2(N[1.7, 2.7], N(1))
    s1 = Singleton(N[1, 2])
    s2 = Singleton(N[2, 2])
    subset, point = ⊆(b1, s1, true)
    @test b1 ⊈ s1 && !subset && point ∈ b1 && point ∉ s1
    @test b2 ⊆ s1 && ⊆(b2, s1, true)[1]
    subset, point = ⊆(b1, s2, true)
    @test b1 ⊈ s2 && !subset && point ∈ b1 && point ∉ s2
    @test b1 ⊈ BallInf(N[1, 2], N(1))
    @test b2 ⊆ BallInf(N[1, 2], N(2))
    subset, point = ⊆(b1, b3, true)
    @test b1 ⊈ b3 && !subset && point ∈ b1 && point ∉ b3
    @test b3 ⊆ b1 && ⊆(b3, b1, true)[1]
    # inclusion of a Ball2 in a polyhedron
    b = Ball2(zeros(N, 4), N(1))
    p = BallInf(zeros(N, 2), N(1)) × BallInf(zeros(N, 2), N(1))
    subset, point = ⊆(b, p, true)
    @test subset && b ⊆ p && point == N[]
    p = BallInf(ones(N, 2), N(1)) × BallInf(zeros(N, 2), N(1))
    subset, point = ⊆(b, p, true)
    @test !subset && b ⊈ p && point ∈ b && point ∉ p

    # intersection
    b1 = Ball2(N[0, 0], N(2))
    b2 = Ball2(N[2, 2], N(2))
    b3 = Ball2(N[4, 4], N(2))
    b4 = Ball2(N[1, 1], N(1))
    intersection_empty, point = is_intersection_empty(b1, b2, true)
    @test !is_intersection_empty(b1, b2) && !intersection_empty && point ∈ b1 && point ∈ b2
    @test is_intersection_empty(b1, b3) && is_intersection_empty(b1, b3, true)[1]
    intersection_empty, point = is_intersection_empty(b2, b1, true)
    @test !is_intersection_empty(b2, b1) && !intersection_empty && point ∈ b2 && point ∈ b1
    intersection_empty, point = is_intersection_empty(b1, b4, true)
    @test !is_intersection_empty(b1, b4) && !intersection_empty && point ∈ b1 && point ∈ b4

    # uniform sampling
    B = Ball2(N[1, 2, 3], N(0.5))
    s = sample(B, 100) # get 100 random elements in B
    @test all(si -> si ∈ B, s)
    # sampling without arguments returns a single vector
    s = sample(B)
    @test s isa AbstractVector{N} && s ∈ B

    # Chebyshev center
    c, r = chebyshev_center_radius(B)
    @test c == center(B) && r == B.radius

    # area/volume
    B = Ball2(zeros(N, 2), N(2))
    @test area(B) == volume(B) == N(pi) * radius(B)^2
    B = Ball2(zeros(N, 3), N(2))
    @test_throws AssertionError area(B)
    @test volume(B) == 4 / 3 * N(pi) * radius(B)^3

    # projection
    b4 = Ball2(N[4, 3, 2, 1], N(2))
    @test project(b4, [2, 4]) == Ball2(N[3, 1], N(2))

    # low/high/extrema
    B = Ball2(N[1, 2], N(2))
    l1 = low(B, 1)
    l2 = low(B, 2)
    h1 = high(B, 1)
    h2 = high(B, 2)
    @test extrema(B, 1) == (l1, h1)
    @test l1 == N(-1) && l2 == N(0) && h1 == N(3) && h2 == N(4)
    @test box_approximation(B) ≈ Hyperrectangle(; low=[l1, l2], high=[h1, h2])
    @test low(B) == [N(-1), N(0)]
    @test high(B) == [N(3), N(4)]
    @test extrema(B) == (low(B), high(B))

    # reflect
    @test reflect(b4) == Ball2(N[-4, -3, -2, -1], N(2))

    # scale
    B = Ball2(N[-2, 3], N(1))
    @test scale(N(2), B) == Ball2(N[-4, 6], N(2))
    @test scale(N(-2), B) == Ball2(N[4, -6], N(2))
end
