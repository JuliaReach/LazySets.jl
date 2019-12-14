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

    # boundedness
    @test isbounded(b)

    # isempty
    @test !isempty(b)

    # isuniversal
    answer, w = isuniversal(b, true)
    @test !isuniversal(b) && !answer && w ∉ b

    # an_element & membership function
    @test an_element(b) ∈ b

    # translation
    @test translate(b, N[1, 2]) == Ball1(N[1, 2], N(2))

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
    @test ispermutation(clist, [HalfSpace(N[1, 1], N(1)),  # x + y <= 1
                                HalfSpace(N[-1, 1], N(1)),  # -x + y <= 1
                                HalfSpace(N[1, -1], N(1)), # x - y <= 1
                                HalfSpace(N[-1, -1], N(1))])  # x + y >= -1
end
