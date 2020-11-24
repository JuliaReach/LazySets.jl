for N in [Float64, Rational{Int}, Float32]
    # random hyperrectangle
    rand(Hyperrectangle)

    # constructor with mixed vectors
    Hyperrectangle(sparsevec([1], N[1], 1), N[1])
    Hyperrectangle(N[1], sparsevec([1], N[1], 1))

    # center/radius/high/low/generators
    c = N[0, 0]
    r = N[1, 1]
    h = Hyperrectangle(c, r)
    @test center(h) == c
    @test radius_hyperrectangle(h) == r
    @test high(h) == r
    @test low(h) == -r
    for i in 1:2
        @test radius_hyperrectangle(h, i) == r[i]
        @test high(h, i) == r[i]
        @test low(h, i) == -r[i]
    end
    gens = [Vector(g) for g in generators(h)]
    @test ispermutation(gens, [N[1, 0], N[0, 1]])
    @test genmat(h) ∈ [N[1 0; 0 1], N[0 1; 1 0]]
    @test ngens(h) == 2
    h_flat = Hyperrectangle(N[1, 2, 3, 4, 5], N[1, 0, 2, 0, 3])
    @test collect(generators(h_flat)) ==
        [SingleEntryVector(1, 5, N(1)), SingleEntryVector(3, 5, N(2)),
         SingleEntryVector(5, 5, N(3))]
    @test genmat(h_flat) == N[1 0 0; 0 0 0; 0 2 0; 0 0 0; 0 0 3]
    @test ngens(h_flat) == 3

    # 1D Hyperrectangle
    h = Hyperrectangle(N[0], N[1])
    # Test Dimension
    @test dim(h) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, h) == N[1]
    d = N[-1]
    @test σ(d, h) == N[-1]

    # 2D Hyperrectangle
    h = Hyperrectangle(N[0, 0], N[1, 1])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, h) == N[1, 1]
    d = N[-1, 1]
    @test σ(d, h) == N[-1, 1]
    d = N[-1, -1]
    @test σ(d, h) == N[-1, -1]
    d = N[1, -1]
    @test σ(d, h) == N[1, -1]

    # 2D Hyperrectangle not 0-centered
    h = Hyperrectangle(N[1, 2], N[1, 1])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, h) == N[2, 3]
    d = N[-1, 1]
    @test σ(d, h) == N[0, 3]
    d = N[-1, -1]
    @test σ(d, h) == N[0, 1]
    d = N[0, -1]
    @test σ(d, h) == N[2, 1]

    # 2D Hyperrectangle not same radius in each direction
    h = Hyperrectangle(N[0, 0], N[1, 2])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, h) == N[1, 2]
    d = N[-1, 1]
    @test σ(d, h) == N[-1, 2]
    d = N[-1, -1]
    @test σ(d, h) == N[-1, -2]
    d = N[1, -1]
    @test σ(d, h) == N[1, -2]

    # 2D Hyperrectangle not centered, not same radius, for vertex representation,
    # radius, and diameter
    h = Hyperrectangle(N[3, 2], N[2, 1])
    vl = vertices_list(h)
    # Test Vertices
    @test ispermutation(vl, [N[1, 1], N[1, 3], N[5, 1], N[5, 3]])
    # norm
    @test norm(h) == norm(N[5, 3], Inf)
    # radius
    @test radius(h) == norm(N[2, 1], Inf)
    # diameter
    @test diameter(h) == norm(N[5, 3] - N[1, 1], Inf)

    # alternative constructor
    c = ones(N, 2)
    r = N[2, 3]
    l = N[-1, -2]
    h = N[3, 4]
    H1 = Hyperrectangle(c, r)
    H2 = Hyperrectangle(low=l, high=h)
    @test H1.center == H2.center
    @test H1.radius == H2.radius
    @test_throws AssertionError Hyperrectangle(low=h, high=l)

    # constructor without bounds check
    @test_throws AssertionError Hyperrectangle(low=N[1], high=N[0]) # default: true
    @test_throws AssertionError Hyperrectangle(low=N[1], high=N[0], check_bounds=true)
    @test Hyperrectangle(low=N[1], high=N[0], check_bounds=false) isa Hyperrectangle

    # Test low and high methods for a hyperrectangle
    H = Hyperrectangle(low=l, high=h)
    @test low(H) == l
    @test high(H) == h

    # boundedness
    @test isbounded(H)

    # isempty
    @test !isempty(H)

    # isuniversal
    answer, w = isuniversal(H, true)
    @test !isuniversal(H) && !answer && w ∉ H

    # membership
    H = Hyperrectangle(N[1, 1], N[2, 3])
    @test N[-1.1, 4.1] ∉ H
    @test N[-1, 4] ∈ H

    # "robust" membership of the support vector
    c = N[1.68, 0.73, 0.64]
    r = N[0.46, 0.24, 1.38]
    H = Hyperrectangle(c, r)
    @test σ(ones(N, 3), H) ∈ H

    # an_element function
    H = Hyperrectangle(N[1, 2], N[3, 4])
    @test an_element(H) ∈ H

    # isflat
    @test isflat(Hyperrectangle(N[1, 2], N[3, 0])) &&
          !isflat(Hyperrectangle(N[1, 2], N[3, 4]))

    # split
    @test ispermutation(split(H, [2, 2]),
        [Hyperrectangle(N[-0.5, 0], N[1.5, 2]),
         Hyperrectangle(N[2.5, 0], N[1.5, 2]),
         Hyperrectangle(N[-0.5, 4], N[1.5, 2]),
         Hyperrectangle(N[2.5, 4], N[1.5, 2])])

    # subset
    H1 = Hyperrectangle(N[1, 3], N[0.5, 0.5])
    H2 = Hyperrectangle(N[2, 2.5], N[0.5, 0.5])
    H3 = Hyperrectangle(N[2, 2], N[2, 3])
    B1 = BallInf(N[2, 2.5], N(0.5))
    B2 = BallInf(N[2, 2], N(1))
    @test H1 ⊈ H2 && H1 ⊆ H3 && H2 ⊆ H3
    subset, point = ⊆(H1, H2, true)
    @test !subset && point ∈ H1 && point ∉ H2
    subset, point = ⊆(H2, H1, true)
    @test !subset && point ∈ H2 && point ∉ H1
    subset, point = ⊆(H1, H3, true)
    @test subset
    @test H2 ⊆ B1 && B1 ⊆ H2
    @test B1 ⊆ B2 && B2 ⊈ B1

    # intersection & intersection emptiness
    H1 = Hyperrectangle(N[1, 1], N[2, 2])
    H2 = Hyperrectangle(N[3, 3], N[2, 2])
    B1 = BallInf(N[2, 4], N(0.5))
    B2 = BallInf(N[3, 3], N(0.5))
    for (X1, X2) in [(H1, H2), (H2, H1)]
        intersection_empty, point = is_intersection_empty(X1, X2, true)
        cap = intersection(X1, X2)
        @test cap isa Hyperrectangle{N} && center(cap) == N[2, 2] &&
              radius_hyperrectangle(cap) == N[1, 1]
        @test !is_intersection_empty(X1, X2) &&
              !intersection_empty && point ∈ X1 && point ∈ X2
    end
    cap = intersection(H1, B1)
    @test cap isa EmptySet{N}
    @test is_intersection_empty(H1, B1) && is_intersection_empty(H1, B1, true)[1]
    cap = intersection(H1, B2)
    @test cap == Hyperrectangle(N[2.75, 2.75], N[0.25, 0.25])
    intersection_empty, point = is_intersection_empty(H1, B2, true)
    @test !is_intersection_empty(H1, B2) && !intersection_empty && point ∈ H1 &&
          point ∈ B2

    # linear map (concrete)
    P = linear_map(N[1 0; 0 2], H1)
    @test P isa Zonotope
    P = linear_map(Diagonal(N[1, 2, 3, 4]), overapproximate(H1 * H1))
    @test P isa Zonotope

    # check that vertices_list is computed correctly if the hyperrectangle
    # is "degenerate"/flat, i.e., its radius contains zeros
    # these tests would crash if all 2^100 vertices were computed (see #92)
    H = Hyperrectangle(ones(N, 100), zeros(N, 100))
    vl = vertices_list(H)
    @test vl == [H.center]
    r = zeros(N, 100); r[1] = N(1)
    H = Hyperrectangle(fill(N(1), 100), r)
    vl = vertices_list(H)
    @test ispermutation(vl, [H.center + r, H.center - r])

    # transform hyperrectangle into a polygon
    H1pol = convert(HPolygon, H1)
    vlist = vertices_list(H1pol)
    @test ispermutation(vlist, [N[3, 3], N[3, -1], N[-1, -1], N[-1, 3]])

    # test that we can produce the list of constraints
    clist = constraints_list(H1)
    @test ispermutation(clist,
                        [HalfSpace(N[1, 0], N(3)), HalfSpace(N[0, 1], N(3)),
                         HalfSpace(N[-1, 0], N(1)), HalfSpace(N[0, -1], N(1))])

    # conversion to and from IntervalArithmetic's IntervalBox type
    B = IntervalBox(IntervalArithmetic.Interval(0, 1), IntervalArithmetic.Interval(0, 1))
    H = convert(Hyperrectangle, B)
    @test convert(IntervalBox, H) == B

    # conversion from other hyperrectangular sets
    H = Hyperrectangle(N[1], N[1])
    @test convert(Hyperrectangle, BallInf(N[1], N(1))) == H
    if VERSION >= v"1.1" || N == Float64
        # in pre-v1.1 versions, IntervalArithmetic always returned a Float64
        # vector in the `center(·)` implementation, leading to an error
        @test convert(Hyperrectangle, Interval(N(0), N(2))) == H
    end

    @test convert(Hyperrectangle, SymmetricIntervalHull(Singleton(N[1, 1]))) ==
          Hyperrectangle(N[0, 0], N[1, 1])

    # conversion to a zonotope
    H = Hyperrectangle(low=N[5//10, 6//10], high=N[124//10, 2355//10])
    Hz = convert(Zonotope, H)
    @test Hz == Zonotope(N[129//20, 2361//20], N[119//20 0//1; 0//1 2349//20])

    # conversion of a hyperrectangle with static array components to a zonotope
    # the specialized method for 2D static arrays is also tested
    H = Hyperrectangle(low=SA[N(5//10), N(6//10)], high=N[N(124//10), N(2355//10)])
    Hz = convert(Zonotope, H)
    Hz_sp = LazySets._convert_2D_static(Zonotope, H)
    Z = Zonotope(SA[N(129//20), N(2361//20)], SA[N(119//20) N(0//1); N(0//1) N(2349//20)])
    @test Hz == Z && Hz_sp == Z

    # rectification
    H = Hyperrectangle(N[-1, 2], N[4, 5])
    Hrect = rectify(H)
    @test Hrect.center == N[1.5, 3.5] &&  Hrect.radius == [1.5, 3.5]

    # Minkowski sum
    H1 = Hyperrectangle(N[0, 1], N[2, 3])
    H2 = Hyperrectangle(N[3, 2], N[1, 0])
    @test minkowski_sum(H1, H2) == Hyperrectangle(N[3, 3], N[3, 3])

    # set difference
    h = Hyperrectangle(low=N[0], high=N[1])
    q = Hyperrectangle(low=N[0], high=N[0.5])
    @test convert(Interval, difference(h, q).array[1]) == Interval(N(0.5), N(1))

    # another set difference test in higher-dimensions
    b = BallInf(N[0, 0, 0], N(1))
    @test isempty(difference(b, b))

    # volume of a hyperrectangular set
    b = BallInf(N[0, 0, 0], N(1))
    @test volume(b) ≈ N(8)

    # concete projection
    @test project(BallInf(zeros(N, 4), N(1)), [1, 2]) == Hyperrectangle(zeros(N, 2), ones(N, 2))
end
