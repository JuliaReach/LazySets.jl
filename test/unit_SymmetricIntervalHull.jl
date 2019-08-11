for N in [Float64, Rational{Int}, Float32]
    # constructor from empty set
    E = EmptySet{N}()
    @test SymmetricIntervalHull(E) == E

    # normal constructor
    h = SymmetricIntervalHull(Ball1(N[2, 3], N(4)))

    # constructor attempt from an unbounded set
    @test_throws AssertionError SymmetricIntervalHull(HalfSpace(N[1], N(0)))

    # dimension
    @test dim(h) == 2

    # support vector
    d = N[1, 0]
    @test σ(d, h)[1] == N(6) && σ(-d, h)[1] == N(-6)
    d = N[0, 1]
    @test σ(d, h)[2] == N(7) && σ(-d, h)[2] == N(-7)

    # boundedness
    @test isbounded(h)

    # isempty
    @test !isempty(h)

    # an_element function
    @test an_element(h) == N[0, 0]

    # radius_hyperrectangle function (uncached and cached)
    h = SymmetricIntervalHull(Ball1(N[2, 3], N(4)))
    for i in 1:2
        @test radius_hyperrectangle(h, 1) == 6
        @test radius_hyperrectangle(h) == N[6, 7]
    end

    # high and low
    h = SymmetricIntervalHull(Singleton(N[1, 2]))
    @test high(h) == N[1, 2]
    @test low(h) == N[-1, -2]
end
