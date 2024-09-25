# Note: in pre-v1.1 versions, IntervalArithmetic always operated with Float64
for N in [Float64, Float32, Rational{Int}]
    # random interval
    rand(Interval)

    # constructor from IntervalArithmetic.Interval
    x = Interval(IA.interval(N(0), N(1)))

    # constructor from a vector
    x = Interval(N[0, 1])

    # constructor from a number
    x = Interval(N(0))
    @test x == Interval(N[0, 0])

    # type-less constructor
    x = Interval(N(0), N(1))

    # constructor with promotion
    y = Interval(0, N(1))
    @test y == x

    # unbounded intervals are caught
    @test_throws AssertionError Interval(-Inf, zero(N))
    @test_throws AssertionError Interval(zero(N), Inf)

    @test dim(x) == 1
    @test center(x) == N[0.5]
    @test center(x, 1) == N(0.5)
    @test_throws AssertionError center(x, 2)
    @test min(x) == N(0) && max(x) == N(1)
    v = vertices_list(x)
    @test N[0] in v && N[1] in v

    # vertices list for degenerate interval
    @test vertices_list(Interval(N(0), N(0))) == [[N(0)]]

    # test interface method an_element and membership
    @test an_element(x) ∈ x
    # number in interval is invalid
    @test_throws MethodError N(1) ∈ x
    # test containment
    @test x ⊆ x
    @test x ⊈ N(0.2) * x
    @test x ⊆ N(2) * x
    @test x ⊆ Interval(N(0), N(2))
    @test x ⊈ Interval(N(-1), N(0.5))

    # concrete linear map
    @test linear_map(hcat(N(2)...), x) == Interval(N(2) * x.dat)
    # alias for scale
    @test linear_map(N(2), x) == Interval(N(2) * x.dat)

    # concrete linear map with zonotope output
    M2 = hcat(N[1, 2, 3])
    @test linear_map(M2, x) == Zonotope(N[0.5, 1.0, 1.5], [N[0.5, 1.0, 1.5]])

    # concrete scale of interval
    @test scale(N(0.5), x) == Interval(N(0.5) * min(x), N(0.5) * max(x))

    # radius_hyperrectangle
    @test radius_hyperrectangle(x) == [N(0.5)]
    @test radius_hyperrectangle(x, 1) == N(0.5)

    # + operator = lazy Minkowski sum of intervals
    y = Interval(N(-2), N(0.5))
    m = x + y
    @test m isa MinkowskiSum
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]
    m = convert(Interval, m) # concretize into an interval
    @test min(m) == N(-2) && max(m) == N(1.5)
    v = vertices_list(m)
    @test N[1.5] in v && N[-2] in v

    # the concrete Minkowski sum of intervals returns an interval
    @test minkowski_sum(x, y) == Interval(N(-2), N(1.5))

    # ispolyhedral
    @test ispolyhedral(x)

    # isempty
    @test !isempty(x)

    # isuniversal
    answer, w = isuniversal(x, true)
    @test !isuniversal(x) && !answer && w ∉ x

    # translation
    @test translate(x, N[2]) == Interval(N(2), N(3))

    # Minkowski sum (test that we get the same results as the concrete operation)
    m = x ⊕ y
    @test m isa MinkowskiSum
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]

    # boundedness
    @test isbounded(x)

    # cartesian product
    cp = x × y
    @test cp isa CartesianProduct
    @test dim(cp) == 2

    # conversion to hyperrectangle
    h = convert(Hyperrectangle, x)
    @test h isa Hyperrectangle && center(h) == radius_hyperrectangle(h) == N[0.5]

    # diameter
    x = Interval(N(1), N(3))
    @test diameter(x) == diameter(x, Inf) == diameter(x, 2) == N(2)

    # split
    intervals = [Interval(N(1), N(3 // 2)), Interval(N(3 // 2), N(2)),
                 Interval(N(2), N(5 // 2)), Interval(N(5 // 2), N(3))]
    @test split(x, 4) == split(x, [4]) == intervals
    @test_throws AssertionError split(x, 0)
    @test_throws AssertionError split(x, [4, 4])

    # concrete intersection
    A = Interval(N(5), N(7))
    B = Interval(N(3), N(6))
    C = intersection(A, B)
    @test C isa Interval
    @test min(C) == N(5) && max(C) == N(6)
    # check empty intersection
    E = intersection(A, Interval(N(0), N(1)))
    @test isempty(E)
    # intersection with half-space
    i = Interval(N(1), N(2))
    hs = HalfSpace(N[1], N(1.5))
    @test intersection(i, hs) == Interval(N(1), N(1.5))
    hs = HalfSpace(N[-2], N(-5))
    @test intersection(i, hs) == EmptySet{N}(1)
    hs = HalfSpace(N[2], N(5))
    @test intersection(i, hs) == i
    # intersection with hyperplane
    hp = Hyperplane(N[2], N(3))
    @test intersection(i, hp) == Singleton(N[1.5])
    hp = Hyperplane(N[-1], N(-3))
    @test intersection(i, hp) == EmptySet{N}(1)
    # other intersections
    Y = Ball1(N[2], N(0.5))
    @test intersection(i, Y) == Interval(N(1.5), N(2))
    Y = ConvexHull(Singleton(N[-5]), Singleton(N[-1]))
    @test intersection(i, Y) == EmptySet{N}(1)

    # disjointness check
    @test !isdisjoint(A, B)
    @test isdisjoint(A, Interval(N(-1), N(1)), false)
    s, w = isdisjoint(A, B, true)
    @test s == false && w ∈ A && w ∈ B
    # disjoint with tolerance
    if N == Float64
        @test !isdisjoint(Interval(0.0, 1.0), Interval(1.0 + 1e-10, 2.0))
        LazySets.set_rtol(Float64, 1e-10)
        @test isdisjoint(Interval(0.0, 1.0), Interval(1.0 + 1e-10, 2.0))
        LazySets.set_rtol(Float64, LazySets.default_tolerance(Float64).rtol) # restore
    end

    # conversion from a hyperrectangular set to an interval
    H = Hyperrectangle(N[0], N[1 / 2])
    A = convert(Interval, H)
    @test A isa Interval && low(A) == [N(-1 / 2)] && high(A) == [N(1 / 2)]

    # conversion from a LazySet to an interval
    M = hcat(N[2])
    B = convert(Interval, M * H)
    @test B isa Interval && low(B) == [N(-1)] && high(B) == [N(1)]
    # conversion to and from IntervalArithmetic.Interval
    B2 = convert(IA.Interval, M * H)
    @test B2 == B.dat
    B3 = convert(Interval, B2)
    @test B3 == B
    # conversion from a non-convex set fails
    U = UnionSet(Interval(1, 2), Interval(3, 4))
    @test_throws AssertionError convert(Interval, U)
    @test_throws AssertionError convert(IA.Interval, U)

    # set difference
    A = Interval(N(5), N(8))
    B = Interval(N(6), N(8))
    C = Interval(N(9), N(10))
    D = Interval(N(6), N(7))
    dAB = difference(A, B)
    dAC = difference(A, C)
    dAD = difference(A, D)
    @test dAB == Interval(N(5), N(6))
    @test dAC == Interval(N(5), N(8))
    @test dAD == UnionSet(Interval(N(5), N(6)), Interval(N(7), N(8)))

    # check if an interval is flat, i.e. if its endpoints coincide (to numerical precision)
    ztol = LazySets._ztol(N) # pick up default absolute zero tolerance value
    @test isflat(Interval(N(0), ztol))
    if N <: AbstractFloat
        @test !isflat(Interval(N(0), 2 * ztol))
    elseif N == Rational{Int}
        @test isflat(Interval(N(0), 2 * ztol))
    end

    # rectification
    x = Interval(N(-2), N(-1))
    @test rectify(x) == Interval(N(0), N(0))
    x = Interval(N(-2), N(2))
    @test rectify(x) == Interval(N(0), N(2))
    x = Interval(N(1), N(2))
    @test rectify(x) == x

    # list of vertices of IA types
    b = IntervalBox(IA.interval(0, 1), IA.interval(0, 1))
    vlistIB = vertices_list(b)
    @test is_cyclic_permutation(vlistIB,
                                [SA[N(1), N(1)], SA[N(0), N(1)], SA[N(1), N(0)], SA[N(0), N(0)]])

    vlistI = vertices_list(b[1])
    @test is_cyclic_permutation(vlistI, [SA[N(0)], SA[N(1)]])

    # generators
    @test ngens(x) == 1
    @test collect(generators(x)) == [N[1 / 2]]
    @test genmat(x) == hcat(N[1 / 2])
    x_degenerate = Interval(N(1), N(1))
    @test ngens(x_degenerate) == 0
    gens = genmat(x_degenerate)
    @test gens == Matrix{N}(undef, 1, 0) && gens isa Matrix{N}
    gens = collect(generators(x_degenerate))
    @test isempty(gens) && gens isa Vector{SingleEntryVector{N}}

    # is_interior_point
    v1 = N[3 / 2]
    v2 = N[1]
    if N <: AbstractFloat
        @test is_interior_point(v1, x) && !is_interior_point(v2, x)
    else
        @test_throws AssertionError is_interior_point(v1, x)
        @test is_interior_point(v1, x; ε=1 // 100) && !is_interior_point(v2, x; ε=1 // 100)
    end
    # different numeric type
    v1 = Float16[3 / 2]
    if N <: AbstractFloat
        v2 = Float16[1]
        @test is_interior_point(v1, x) && !is_interior_point(v2, x)
    else
        @test_throws ArgumentError is_interior_point(v1, x)
        @test_throws ArgumentError is_interior_point(v1, x; ε=1 // 100)
    end

    # Chebyshev center
    c, r = chebyshev_center_radius(x)
    @test c == center(x) && r == N(1 // 2)

    # reflect
    @test reflect(Interval(N(1), N(2))) == Interval(N(-2), N(-1))
    @test reflect(Interval(N(-1), N(2))) == Interval(N(-2), N(1))
    @test reflect(Interval(N(-2), N(-1))) == Interval(N(1), N(2))

    # issubset
    I13 = Interval(N(1), N(3))
    I02 = Interval(N(0), N(2))
    I24 = Interval(N(2), N(4))
    I03 = Interval(N(0), N(3))
    I14 = Interval(N(1), N(4))
    for I in (I02, I24)
        res, w = ⊆(I13, I, true)
        @test !(⊆(I13, I)) && !res && w ∈ I13 && w ∉ I
    end
    res, w = ⊆(I13, I13, true)
    @test ⊆(I13, I13) && res && w == N[]
    # isstrictsubset
    for I in (I02, I24, I13)
        res, w = ⊂(I13, I, true)
        @test !(⊂(I13, I)) && !res && w == N[]
    end
    for I in (I03, I14)
        res, w = ⊂(I13, I, true)
        @test ⊂(I13, I) && res && w ∈ I && w ∉ I13
    end

    # permute
    @test permute(x, [1]) == x
    @test_throws AssertionError permute(x, Int[])
    @test_throws AssertionError permute(x, Int[2])
    @test_throws AssertionError permute(x, Int[1, 1])

    # project
    @test project(x, [1]) == x
    @test_throws AssertionError project(x, Int[])
    @test_throws AssertionError project(x, Int[2])
    @test_throws AssertionError project(x, Int[1, 1])

    # isapprox
    @test x ≈ x ≈ translate(x, [1e-8])
    @test !(x ≈ translate(x, [1e-4]))

    # convex_hull & linear_combination
    I1 = Interval(N(0), N(1))
    I2 = Interval(N(2), N(3))
    @test convex_hull(I1, I2) == linear_combination(I1, I2) == Interval(N(0), N(3))

    # complement
    C = complement(I2)
    L = HalfSpace(N[1], N(2))
    H = HalfSpace(N[-1], N(-3))
    @test length(C) == 2 && ispermutation([C[1], C[2]], [L, H])

    # low, high, extrema
    @test extrema(I1) == (low(I1), high(I1)) == (N[0], N[1])
    @test extrema(I1, 1) == (low(I1, 1), high(I1, 1)) == (N(0), N(1))
    @test_throws AssertionError low(I1, 2)
    @test_throws AssertionError high(I1, 2)
    @test_throws AssertionError extrema(I1, 2)

    # norm
    @test norm(I1) == N(1)
    @test norm(Interval(N(-2), N(1))) == N(2)

    # volume
    @test volume(I1) == 1
    @test volume(Interval(N(-2), N(1))) == N(3)

    # affine_map
    M = hcat(N[2])
    v = N[-3]
    @test affine_map(M, I1, v) == Interval(N(-3), N(-1))
    M2 = N[2; -2;;]
    v = N[-1, -1]
    H = affine_map(M2, I1, v)
    @test H isa Zonotope && isequivalent(H, LineSegment(N[-1, -1], N[1, -3]))

    # exponential_map
    @test exponential_map(M, I1) == Interval(N(0), N(exp(N(2))))

    # isequivalent
    @test isequivalent(I1, I1) && !isequivalent(I1, I2)

    # distance
    @test distance(I1, I2) == distance(I2, I1) == N(1)
    I3 = Interval(N(1 // 2), N(2))
    @test distance(I1, I3) == distance(I2, I3) == distance(I1, I1) == N(0)
end
