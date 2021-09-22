for N in [Float64, Float32]
    # random ball
    rand(Ballp)

    # 1D Ball3
    b = Ballp(N(3), N[0], N(1))
    # dimension
    @test dim(b) == 1
    # support vector
    d = N[1]
    @test σ(d, b) == N[1]
    d = N[-1]
    @test σ(d, b) == N[-1]

    # 2D Ball3
    b = Ballp(N(3), N[0, 0], N(2))
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
    d = N[0, 0]
    @test σ(d, b) == N[0, 0]

    # constructors that fall back to specialized set types
    @test Ballp(N(1), N[3], N(4)) isa Ball1
    @test Ballp(N(2), N[3], N(4)) isa Ball2
    @test Ballp(N(Inf), N[3], N(4)) isa BallInf

    # center
    @test center(b) == b.center
    @test an_element(b) isa AbstractVector{N}

    # boundedness
    @test isbounded(b)

    # isempty
    @test !isempty(b)

    # isuniversal
    answer, w = isuniversal(b, true)
    @test !isuniversal(b) && !answer && w ∉ b

    # membership & an_element
    @test an_element(b) ∈ b

    # translation
    @test translate(b, N[1, 2]) == Ballp(N(3), N[1, 2], N(2))
    bb = Ballp(N(3), N[0, 0], N(1))
    @test translate!(bb, N[1, 1]) == Ballp(N(3), N[1, 1], N(1)) == bb

    # projection
    b4 = Ballp(N(3), N[4, 3, 2, 1], N(2))
    @test project(b4, [2, 4]) == Ballp(N(3), N[3, 1], N(2))
end
