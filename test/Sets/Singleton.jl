using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation
using LazySets.ReachabilityBase.Arrays: SingleEntryVector
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # constructors
    S1 = Singleton(N[1, 2])
    S2 = Singleton(N(1), N(2))
    @test S1 == S2

    # isoperationtype
    @test !isoperationtype(Singleton)

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

    # ispolyhedral
    @test ispolyhedral(s)

    # element function
    @test element(s) == s.element
    for i in 1:2
        @test element(s, i) == s.element[i]
    end

    # isempty
    @test !isempty(s)

    # isuniversal
    answer, w = isuniversal(s, true)
    @test !isuniversal(s) && !answer && w ∉ s

    # sampling
    @test sample(s, 2) == [element(s), element(s)]

    # membership
    S = Singleton(N[1, 1])
    @test N[0.9, 1.1] ∉ S
    @test N[1, 1] ∈ S

    # check membership error message
    @test_throws ArgumentError S.element ⊆ S

    # an_element function
    @test an_element(S) == element(S)

    # vertices / vertices_list
    @test collect(vertices(S)) == vertices_list(S)
    @test vertices_list(S) == [element(S)]

    # radius_hyperrectangle
    @test iszero(radius_hyperrectangle(S))

    # high and low
    @test high(S) == element(S)
    @test low(S) == element(S)

    # constraints_list
    clist4 = constraints_list(S)
    @test ispermutation(clist4,
                        [HalfSpace(N[1, 0], N(1)), HalfSpace(N[0, 1], N(1)),
                         HalfSpace(N[-1, 0], N(-1)), HalfSpace(N[0, -1], N(-1))])
    clist3 = constraints_list(S; min_constraints=true)
    @test length(clist3) == 3 && isequivalent(HPolygon(clist4), HPolygon(clist3))

    # concrete linear map
    M = N[0 1; -1 0]
    @test element(linear_map(M, S)) == an_element(M * S)

    # translation
    @test translate(s, N[1, 2]) == Singleton(N[2, 4])
    ss = copy(s)
    @test translate!(ss, N[1, 2]) == Singleton(N[2, 4]) == ss

    # subset
    s1 = Singleton(N[0, 1])
    s2 = Singleton(N[0, 3])
    p1 = VPolygon([N[0, 0], N[0, 2]])
    p2 = VPolygon([N[0, 0], N[0, 2], N[2, 0]])
    b = BallInf(N[0, 1], N(1))
    @test s1 ⊆ p1 && ⊆(s1, p1, true)[1]
    subset, point = ⊆(s2, p2, true)
    @test s2 ⊈ p2 && !subset && point ∈ s2 && point ∉ p2
    @test s1 ⊆ s1 && ⊆(s1, s1, true)[1]
    subset, point = ⊆(s1, s2, true)
    @test s1 ⊈ s2 && !subset && point ∈ s1 && point ∉ s2
    subset, point = ⊆(s1, b, true)
    @test subset && s1 ⊆ b
    subset, point = ⊆(s2, b, true)
    @test s2 ⊈ b && !subset && point ∈ s2 && point ∉ b

    # intersection
    S1 = Singleton(N[1, 1])
    S2 = Singleton(N[0, 0])
    S3 = ZeroSet{N}(2)
    H = BallInf(N[1, 1], N(0.5))
    M = LinearMap(N[1 0; 0 1], H)
    @test intersection(S1, H) == intersection(H, S1) == S1
    @test intersection(S2, H) == intersection(H, S2) == EmptySet{N}(2)
    @test intersection(S3, H) == intersection(H, S3) == EmptySet{N}(2)
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

    # rectification
    @test rectify(S1) == S1
    @test rectify(Singleton(N[-1, -1])) == S2

    # concrete minkowski sum
    @test minkowski_sum(Singleton(N[1, 2]), Singleton(N[3, 4])) == Singleton(N[4, 6])

    # projection
    S = Singleton(N[1, 2, 3])
    @test project(S, 1:2) == Singleton(N[1, 2])
    @test project(S, [1, 3]) == Singleton(N[1, 3])

    # concrete cartesian product
    S1 = Singleton(N[1, 2, 3])
    S2 = Singleton(N[4, 5, 6])
    @test cartesian_product(S1, S2) == Singleton(N[1, 2, 3, 4, 5, 6])
    @test concretize(S1 × S2 × S1) == Singleton(N[1, 2, 3, 4, 5, 6, 1, 2, 3])

    # permutation
    S = Singleton(N[3, 4, 5])
    @test permute(S, 1:3) == S
    @test permute(S, [2, 1, 3]) == Singleton(N[4, 3, 5])

    # Chebyshev center
    c, r = chebyshev_center_radius(S)
    @test c == element(S) && r == zero(N)

    # reflect
    S = Singleton(N[0, -1, 2])
    @test reflect(S) == Singleton(N[0, 1, -2])

    # scale/scale!
    S2 = copy(S)
    scale!(N(2), S2)
    @test scale(N(2), S) == S2 == Singleton(N[0, -2, 4])

    # singleton_list
    @test singleton_list(S) == [S]
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(Singleton; N=N) isa Singleton{N}
end
