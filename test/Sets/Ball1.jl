for N in [Float64, Rational{Int}, Float32]
    # random ball
    rand(Ball1)

    # 1D Ball1
    b = Ball1(N[0], N(1))
    # dimension
    @test dim(b) == 1
    # support vector
    d = N[1]
    @test σ(d, b) == N[1]
    d = N[-1]
    @test σ(d, b) == N[-1]

    # 2D Ball1
    b = Ball1(N[0, 0], N(1))
    # dimension
    @test dim(b) == 2
    # support vector
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

    # 2D Ball1 not 0-centered
    b = Ball1(N[1, 2], N(1))
    # dimension
    @test dim(b) == 2
    # support vector
    d = N[1, 0]
    @test σ(d, b) == N[2, 2]
    d = N[-1, 0]
    @test σ(d, b) == N[0, 2]
    d = N[0, 1]
    @test σ(d, b) == N[1, 3]
    d = N[0, -1]
    @test σ(d, b) == N[1, 1]

    # 2D Ball1 radius != 1
    b = Ball1(N[0, 0], N(2))
    # dimension
    @test dim(b) == 2
    # support vector
    d = N[1, 0]
    @test σ(d, b) == N[2, 0]
    d = N[-1, 0]
    @test σ(d, b) == N[-2, 0]
    d = N[0, 1]
    @test σ(d, b) == N[0, 2]
    d = N[0, -1]
    @test σ(d, b) == N[0, -2]

    # center
    @test center(b) == N[0, 0]

    # radius_ball
    @test LazySets.radius_ball(b) == N(2)

    # ball_norm
    @test LazySets.ball_norm(b) == N(1)

    # low/high/extrema
    @test extrema(b) == (low(b), high(b)) == (N[-2, -2], N[2, 2])
    @test extrema(b, 1) == (low(b, 1), high(b, 1)) == (N(-2), N(2))

    # boundedness
    @test isbounded(b)

    # isoperationtype
    @test !isoperationtype(typeof(b))

    # ispolyhedral
    @test ispolyhedral(b)

    # isempty
    @test !isempty(b)

    # isuniversal
    answer, w = isuniversal(b, true)
    @test !isuniversal(b) && !answer && w ∉ b

    # an_element & membership function
    e = an_element(b)
    answer = e ∈ b
    @test answer == ∈(e, b, true) == ∈(e, b, false) == true
    e += N[2 * b.radius, 1]
    answer = e ∈ b
    @test answer == ∈(e, b, true) == ∈(e, b, false) == false

    # translation
    @test translate(b, N[1, 2]) == Ball1(N[1, 2], N(2))
    bb = Ball1(N[0, 0], N(1))
    @test translate!(bb, N[1, 1]) == Ball1(N[1, 1], N(1)) == bb

    # vertices_list
    vl = vertices_list(b)
    @test ispermutation(vl, [N[2, 0], N[0, 2], N[-2, 0], N[0, -2]])

    # check that vertices_list for zero radius doesn't repeat vertices
    b = Ball1(N[1, 2], N(0))
    vl = vertices_list(b)
    @test vl == [center(b)]

    # list of constraints
    b = Ball1(N[0, 0], N(1))
    clist = constraints_list(b)
    @test ispermutation(clist,
                        [HalfSpace(N[1, 1], N(1)),  # x + y <= 1
                         HalfSpace(N[-1, 1], N(1)),  # -x + y <= 1
                         HalfSpace(N[1, -1], N(1)), # x - y <= 1
                         HalfSpace(N[-1, -1], N(1))])  # x + y >= -1

    # projection
    b4 = Ball1(N[4, 3, 2, 1], N(2))
    @test project(b4, [2, 4]) == Ball1(N[3, 1], N(2))

    # reflect
    @test reflect(b4) == Ball1(N[-4, -3, -2, -1], N(2))

    # scale
    B = Ball1(N[-2, 3], N(1))
    @test scale(N(2), B) == Ball1(N[-4, 6], N(2))
    @test scale(N(-2), B) == Ball1(N[4, -6], N(2))
end
