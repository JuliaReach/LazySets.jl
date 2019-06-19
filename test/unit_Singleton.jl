for N in [Float64, Rational{Int}, Float32]
    # random singleton
    rand(Singleton)

    # center/radius/high/low/generators
    c = N[0, 0]
    r = N[0, 0]
    s = Singleton(c)
    @test center(s) == c
    @test radius_hyperrectangle(s) == r
    @test high(s) == r
    @test low(s) == -r
    for i in 1:2
        @test radius_hyperrectangle(s, i) == r[i]
        @test high(s, i) == r[i]
        @test low(s, i) == -r[i]
    end
    @test collect(generators(s)) == Vector{Vector{N}}()
    @test genmat(s) == Matrix{N}(undef, 2, 0)
    @test ngens(s) == 0

    # 1D singleton
    s = Singleton(N[1])
    d = N[1]
    @test σ(d, s) == N[1]
    d = N[-1]
    @test σ(d, s) == N[1]
    d = N[0]
    @test σ(d, s) == N[1]

    # 2D singleton
    s = Singleton(N[1, 2])
    d = N[1, 0]
    @test σ(d, s) == N[1, 2]
    d = N[-1, 0.5]
    @test σ(d, s) == N[1, 2]
    d = N[0, 0]
    @test σ(d, s) == N[1, 2]

    # test that abstract vectors can be used as well
    S = Singleton(SingleEntryVector(1, 1000, 0.5))
    @test element(S) == SingleEntryVector(1, 1000, 0.5)

    # boundedness
    @test isbounded(s)

    # element function
    @test element(s) == s.element
    for i in 1:2
        @test element(s, i) == s.element[i]
    end

    # isempty
    @test !isempty(s)

    # membership
    S = Singleton(N[1, 1])
    !∈(N[0.9, 1.1], S)
    ∈(N[1, 1], S)

    # check membership error message
    @test_throws ErrorException S.element ⊆ S

    # an_element function
    @test an_element(S) ∈ S

    # vertices_list
    @test vertices_list(S)[1] ∈ S

    # radius_hyperrectangle
    @test iszero(radius_hyperrectangle(S))

    # high and low
    @test high(S) == N[1, 1]
    @test low(S) == N[1, 1]

    # concrete linear map
    M = N[0 1; -1 0]
    @test element(linear_map(M, S)) == an_element(M * S)

    # translation
    @test translate(s, N[1, 2]) == Singleton(N[2, 4])

    # subset
    s1 = Singleton(N[0, 1])
    s2 = Singleton(N[0, 3])
    p1 = VPolygon([N[0,0], N[0, 2]])
    p2 = VPolygon([N[0,0], N[0, 2], N[2, 0]])
    b = BallInf(N[0, 1], N(1))
    @test ⊆(s1, p1) && ⊆(s1, p1, true)[1]
    subset, point = ⊆(s2, p2, true)
    @test !⊆(s2, p2) && !subset && point ∈ s2 && !(point ∈ p2)
    @test ⊆(s1, s1) && ⊆(s1, s1, true)[1]
    subset, point = ⊆(s1, s2, true)
    @test !⊆(s1, s2) && !subset && point ∈ s1 && !(point ∈ s2)
    subset, point = ⊆(s1, b, true)
    @test subset && ⊆(s1, b)
    subset, point = ⊆(s2, b, true)
    @test !⊆(s2, b) && !subset && point ∈ s2 && !(point ∈ b)

    # intersection
    S1 = Singleton(N[1, 1])
    S2 = Singleton(N[0, 0])
    S3 = ZeroSet{N}(2)
    H = BallInf(N[1, 1], N(0.5))
    M = LinearMap(N[1 0; 0 1], H)
    @test intersection(S1, H) == intersection(H, S1) == S1
    @test intersection(S2, H) == intersection(H, S2) == EmptySet{N}()
    @test intersection(S3, H) == intersection(H, S3) == EmptySet{N}()
    @test is_intersection_empty(S1, S2) && is_intersection_empty(S1, S2, true)[1]
    intersection_empty, point = is_intersection_empty(S2, S3, true)
    @test !is_intersection_empty(S2, S3) && !intersection_empty &&
        point ∈ S2 && point ∈ S3
    for X in [H, M]
        intersection_empty, point = is_intersection_empty(S1, X, true)
        @test !is_intersection_empty(S1, X) && !intersection_empty &&
            point ∈ S1 && point ∈ X
        @test is_intersection_empty(S2, X) &&
            is_intersection_empty(S2, X, true)[1]
    end
end
