@static if VERSION >= v"1.1"
    Ns = [Float64, Float32, Rational{Int}]
else
    # in pre-v1.1 versions, IntervalArithmetic always operated with Float64
    Ns = [Float64]
end

for N in Ns
    # random interval
    rand(Interval)

    # constructor from IntervalArithmetic.Interval
    x = Interval(IntervalArithmetic.Interval(N(0), N(1)))

    # constructor from a vector
    x = Interval(N[0, 1])

    # type-less constructor
    x = Interval(N(0), N(1))

    @test dim(x) == 1
    @test center(x) == N[0.5]
    @test min(x) == N(0) && max(x) == N(1)
    v = vertices_list(x)
    @test N[0] in v && N[1] in v
    # test interface method an_element and membership
    @test an_element(x) ∈ x
    # test containment
    @test x ⊆ x
    @test x ⊈ N(0.2) * x
    @test x ⊆ N(2) * x
    @test x ⊆ Interval(N(0), N(2))
    @test x ⊈ Interval(N(-1), N(0.5))

    # concrete linear map
    @test linear_map(hcat(N(2)...), x) == Interval(N(2) * x.dat)

    # concrete linear map with zonotope output
    M2 = hcat(N[1, 2, 3])
    @test linear_map(M2, x) == Zonotope(N[0.5, 1.0, 1.5], [N[0.5, 1.0, 1.5]])

    # concrete scale of interval
    @test scale(N(0.5), x) == Interval(N(0.5)*min(x), N(0.5)*max(x))

    # radius_hyperrectangle
    @test radius_hyperrectangle(x) == [N(0.5)]
    @test radius_hyperrectangle(x, 1) == N(0.5)

    # + operator (= concrete Minkowski sum of intervals)
    y = Interval(N(-2), N(0.5))
    m = x + y
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]
    @test min(m) == N(-2) && max(m) == N(1.5)
    v = vertices_list(m)
    @test N[1.5] in v && N[-2] in v

    # subtraction
    d = x - y
    @test dim(d) == 1
    @test σ(N[1], d) == N[3]
    @test σ(N[-1], d) == N[-0.5]
    @test min(d) == N(-0.5) && max(d) == N(3)
    v = vertices_list(d)
    @test N[-0.5] in v && N[3] in v

    # product of intervals: use the * operator
    p = x * y
    @test dim(p) == 1
    @test σ(N[1], p) == N[0.5]
    @test σ(N[-1], p) == N[-2]
    v = vertices_list(p)
    @test N[0.5] in v && N[-2] in v

    # test different arithmetic operations
    r = (x + y) - (d + p)
    @test min(r) == N(-5.5) && max(r) == N(4)

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
    @test intersection(i, hs) == EmptySet{N}()
    hs = HalfSpace(N[2], N(5))
    @test intersection(i, hs) == i
    # intersection with hyperplane
    hp = Hyperplane(N[2], N(3))
    @test intersection(i, hp) == Singleton(N[1.5])
    hp = Hyperplane(N[-1], N(-3))
    @test intersection(i, hp) == EmptySet{N}()
    # other intersections
    Y = Ball1(N[2], N(0.5))
    @test intersection(i, Y) == Interval(N(1.5), N(2))
    Y = ConvexHull(Singleton(N[-5]), Singleton(N[-1]))
    @test intersection(i, Y) == EmptySet{N}()

    # disjointness check
    @test !isdisjoint(A, B)
    @test isdisjoint(A, Interval(N(-1), N(1)), false)
    s, w = isdisjoint(A, B, true)
    @test s == false && w ∈ A && w ∈ B

    # conversion from a hyperrectangular set to an interval
    H = Hyperrectangle(N[0], N[1/2])
    A = convert(Interval, H)
    @test A isa Interval && low(A) == [N(-1/2)] && high(A) == [N(1/2)]

    # conversion from a lazyset to an interval
    M = hcat(N[2])
    B = convert(Interval, M*H)
    @test B isa Interval && low(B) == [N(-1)] && high(B) == [N(1)]

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
        @test !isflat(Interval(N(0), 2*ztol))
    elseif N == Rational{Int}
        @test isflat(Interval(N(0), 2*ztol))
    end
end
