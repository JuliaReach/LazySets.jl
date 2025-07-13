for N in @tN([Float64, Float32, Rational{Int}])
    # constructor from empty set
    E = EmptySet{N}(2)
    @test SymmetricIntervalHull(E) == E

    # normal constructor
    B = Ball1(N[2, 3], N(4))
    h = SymmetricIntervalHull(B)

    # alias
    @test ⊡(B) == h

    # constructor attempt from an unbounded set
    h0 = HalfSpace(N[1], N(0))
    @test_throws AssertionError SymmetricIntervalHull(h0)

    # pass boundedness constructor flag
    S = SymmetricIntervalHull(h0; check_boundedness=false)
    @test S.X == h0

    # dimension
    @test dim(h) == 2

    # support vector
    d = N[1, 0]
    @test σ(d, h)[1] == N(6) && σ(-d, h)[1] == N(-6)
    d = N[0, 1]
    @test σ(d, h)[2] == N(7) && σ(-d, h)[2] == N(-7)

    # boundedness
    @test isbounded(h) && isboundedtype(typeof(h))

    # ispolyhedral
    @test ispolyhedral(h)

    # isempty
    @test !isempty(h)

    # isuniversal
    answer, w = isuniversal(h, true)
    @test !isuniversal(h) && !answer && w ∉ h

    # an_element function
    @test an_element(h) == N[0, 0]

    # radius_hyperrectangle function (uncached and cached)
    h = SymmetricIntervalHull(B)
    for i in 1:2
        @test radius_hyperrectangle(h, 1) == 6
        @test radius_hyperrectangle(h) == N[6, 7]
    end

    # high and low
    h = SymmetricIntervalHull(Singleton(N[1, 2]))
    @test high(h) == N[1, 2]
    @test low(h) == N[-1, -2]

    # center
    @test center(h) == N[0, 0]
    @test center(h, 1) == zero(N)

    # concretize
    B = Ball1(N[0, 0], N(1))
    h = SymmetricIntervalHull(B)
    @test concretize(h) == symmetric_interval_hull(B)
end
