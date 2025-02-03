function isidentical(::Interval, ::Interval)
    return false
end

function isidentical(X1::Interval{N}, X2::Interval{N}) where {N}
    return X1.dat == X2.dat
end

for N in [Float64, Float32, Rational{Int}]
    # default constructor from IntervalArithmetic.Interval
    itv = IA.interval(N(0), N(1))
    X = Interval(itv)
    @test X isa Interval{N} && X.dat == itv

    # constructors from two numbers, from a vector, and with promotion
    for Y in (Interval(N(0), N(1)), Interval(N[0, 1]), Interval(0, N(1)))
        @test isidentical(Y, X)
    end

    # constructor from a number
    Y = Interval(N(0))
    @test Y isa Interval{N} && Y == Interval(N(0), N(0))

    # unbounded intervals are caught
    @test_throws AssertionError Interval(-Inf, zero(N))
    @test_throws AssertionError Interval(zero(N), Inf)

    # convert
    Y = Interval(N(-1), N(1))
    Z = convert(Interval, Rectification(Y))
    @test isidentical(Z, X)
    Z = convert(Interval, Y + Y)
    @test Z isa Interval{N} && Z == Interval(N(-2), N(2))

    # an_element
    x = an_element(X)
    @test x isa Vector{N} && length(x) == 1 && x[1] isa N && N(0) <= x[1] <= N(1)

    # area
    @test_throws AssertionError area(X)

    # chebyshev_center_radius
    c, r = chebyshev_center_radius(X)
    @test c isa Vector{N} && c == [N(1 // 2)]
    @test r isa N && r == N(1 // 2)

    # complement
    Y = complement(X)
    @test Y isa UnionSet{N} && length(Y) == 2
    @test ispermutation(array(Y), [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))])

    # concretize
    X2 = concretize(X)
    @test isidentical(X, X2)

    # constrained_dimensions
    @test constrained_dimensions(X) == 1:1

    # constraints_list
    clist = constraints_list(X)
    @test ispermutation(clist, [HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(0))])

    # constraints
    clist2 = collect(constraints(X))
    @test ispermutation(clist2, clist)

    # convex_hull (unary)
    X2 = convex_hull(X)
    @test isidentical(X, X2)

    # copy
    X2 = copy(X)
    @test isidentical(X, X2)

    # diameter
    for res in (diameter(X), diameter(X, Inf), diameter(X, 2))
        @test res isa N && res == N(1)
    end

    # dim
    @test dim(X) == 1

    # eltype
    @test eltype(X) == N
    @test eltype(typeof(X)) == N

    # extrema
    res = extrema(X)
    @test res isa Tuple{Vector{N}, Vector{N}} && res[1] == N[0] && res[2] == N[1]
    res = extrema(X, 1)
    @test res isa Tuple{N, N} && res[1] == N(0) && res[2] == N(1)
    @test_throws AssertionError extrema(X, 2)

    # high
    res = high(X)
    @test res isa Vector{N} && res == N[1]
    res = high(X, 1)
    @test res isa N && res == N(1)
    @test_throws AssertionError high(X, 2)

    # isbounded
    @test isbounded(X)

    # isboundedtype
    @test isboundedtype(typeof(X))

    # isconvextype
    @test isconvextype(typeof(X))

    # isempty
    @test !isempty(X)

    # isoperation
    @test !isoperation(X)

    # isoperationtype
    @test !isoperationtype(typeof(X))

    # ispolyhedral
    @test ispolyhedral(X)

    # isuniversal
    res, w = isuniversal(X, true)
    @test !isuniversal(X) && !res && w isa Vector{N} && w ∉ X

    # low
    res = low(X)
    @test res isa Vector{N} && res == N[0]
    res = low(X, 1)
    @test res isa N && res == N(0)
    @test_throws AssertionError low(X, 2)

    # norm
    res = norm(X)
    @test res isa N && res == N(1)
    @test norm(Interval(N(-2), N(1))) == N(2)

    # polyhedron
    if test_suite_polyhedra
        Y = polyhedron(X)
        @test Y isa Polyhedra.Interval{N}
        @test ispermutation(Z.points.points, [[N(0)], [N(1)]])
    end




    # random interval
    rand(Interval)

    @test center(X) == N[0.5]
    @test center(X, 1) == N(0.5)
    @test_throws AssertionError center(X, 2)
    @test min(X) == N(0) && max(X) == N(1)
    v = vertices_list(X)
    @test N[0] in v && N[1] in v

    # vertices list for degenerate interval
    @test vertices_list(Interval(N(0), N(0))) == [[N(0)]]

    # number in interval is invalid
    @test_throws MethodError N(1) ∈ X
    # test containment
    @test X ⊆ X
    @test X ⊈ N(0.2) * X
    @test X ⊆ N(2) * X
    @test X ⊆ Interval(N(0), N(2))
    @test X ⊈ Interval(N(-1), N(0.5))

    # concrete linear map
    @test linear_map(hcat(N(2)...), X) == Interval(N(2) * X.dat)

    # concrete linear map with zonotope output
    M2 = hcat(N[1, 2, 3])
    @test linear_map(M2, X) == Zonotope(N[0.5, 1.0, 1.5], [N[0.5, 1.0, 1.5]])

    # concrete scale of interval
    @test scale(N(0.5), X) == Interval(N(0.5) * min(X), N(0.5) * max(X))

    # radius_hyperrectangle
    @test radius_hyperrectangle(X) == [N(0.5)]
    @test radius_hyperrectangle(X, 1) == N(0.5)

    # the concrete Minkowski sum of intervals returns an interval
    @test minkowski_sum(X, Y) == Interval(N(-2), N(1.5))

    # translation
    @test translate(X, N[2]) == Interval(N(2), N(3))

    # Minkowski sum (test that we get the same results as the concrete operation)
    m = X ⊕ Y
    @test m isa MinkowskiSum
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]

    # cartesian product
    cp = X × Y
    @test cp isa CartesianProduct
    @test dim(cp) == 2

    # conversion to hyperrectangle
    h = convert(Hyperrectangle, X)
    @test h isa Hyperrectangle && center(h) == radius_hyperrectangle(h) == N[0.5]

    # split
    intervals = [Interval(N(1), N(3 // 2)), Interval(N(3 // 2), N(2)),
                 Interval(N(2), N(5 // 2)), Interval(N(5 // 2), N(3))]
    @test split(X, 4) == split(X, [4]) == intervals
    @test_throws AssertionError split(X, 0)
    @test_throws AssertionError split(X, [4, 4])

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
    Y = UnionSet(Interval(1, 2), Interval(3, 4))
    @test_throws AssertionError convert(Interval, Y)
    @test_throws AssertionError convert(IA.Interval, Y)

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
    X = Interval(N(-2), N(-1))
    @test rectify(X) == Interval(N(0), N(0))
    X = Interval(N(-2), N(2))
    @test rectify(X) == Interval(N(0), N(2))
    X = Interval(N(1), N(2))
    @test rectify(X) == X

    # list of vertices of IA types
    b = IntervalBox(IA.interval(0, 1), IA.interval(0, 1))
    vlistIB = vertices_list(b)
    @test is_cyclic_permutation(vlistIB,
                                [SA[N(1), N(1)], SA[N(0), N(1)], SA[N(1), N(0)], SA[N(0), N(0)]])

    vlistI = vertices_list(b[1])
    @test is_cyclic_permutation(vlistI, [SA[N(0)], SA[N(1)]])

    # generators
    @test ngens(X) == 1
    @test collect(generators(X)) == [N[1 / 2]]
    @test genmat(X) == hcat(N[1 / 2])
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
        @test is_interior_point(v1, X) && !is_interior_point(v2, X)
    else
        @test_throws AssertionError is_interior_point(v1, X)
        @test is_interior_point(v1, X; ε=1 // 100) && !is_interior_point(v2, X; ε=1 // 100)
    end
    # different numeric type
    v1 = Float16[3 / 2]
    if N <: AbstractFloat
        v2 = Float16[1]
        @test is_interior_point(v1, X) && !is_interior_point(v2, X)
    else
        @test_throws ArgumentError is_interior_point(v1, X)
        @test_throws ArgumentError is_interior_point(v1, X; ε=1 // 100)
    end

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
    @test permute(X, [1]) == X
    @test_throws AssertionError permute(X, Int[])
    @test_throws AssertionError permute(X, Int[2])
    @test_throws AssertionError permute(X, Int[1, 1])

    # project
    @test project(X, [1]) == X
    @test_throws AssertionError project(X, Int[])
    @test_throws AssertionError project(X, Int[2])
    @test_throws AssertionError project(X, Int[1, 1])

    # isapprox
    @test X ≈ X ≈ translate(X, [1e-8])
    @test !(X ≈ translate(X, [1e-4]))

    # convex_hull & linear_combination
    I1 = Interval(N(0), N(1))
    I2 = Interval(N(2), N(3))
    @test convex_hull(I1, I2) == linear_combination(I1, I2) == Interval(N(0), N(3))

    # volume
    @test volume(I1) == 1
    @test volume(Interval(N(-2), N(1))) == N(3)

    # affine_map
    M = hcat(N[2])
    v = N[-3]
    @test affine_map(M, I1, v) == Interval(N(-3), N(-1))
    M2 = hcat(N[2, -2])
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

for N in [Float64, Float32]
    E = EmptySet{N}(2)

    # rationalize
    E2 = rationalize(E)
    @test E2 isa EmptySet{Rational{Int}} && dim(E2) == 2
    @test_throws MethodError rationalize(E2)

    # is_interior_point
    @test_throws AssertionError is_interior_point(N[0], E)
    @test !is_interior_point(N[0, 0], E)
end
