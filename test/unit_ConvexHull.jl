for N in [Float64, Rational{Int}, Float32]
    # ConvexHull of two 2D Ball1
    b1 = Ball1(N[0, 0], N(1))
    b2 = Ball1(N[1, 2], N(1))
    # Test Construction
    ch = ConvexHull(b1, b2)
    @test ch == CH(b1, b2)

    # swap
    ch2 = swap(ch)
    @test ch.X == ch2.Y && ch.Y == ch2.X

    # Test Dimension
    @test dim(ch) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, ch) == N[2, 2]
    d = N[-1, 0]
    @test σ(d, ch) == N[-1, 0]
    d = N[0, 1]
    @test σ(d, ch) == N[1, 3]
    d = N[0, -1]
    @test σ(d, ch) == N[0, -1]

    # boundedness
    @test isbounded(ch)
    @test !isbounded(CH(Singleton(N[1]), HalfSpace(N[1], N(1))))

    # isempty
    @test !isempty(ch)

    # test convex hull of a set of points using the default algorithm
    points = to_N(N, [[0.9, 0.2], [0.4, 0.6], [0.2, 0.1], [0.1, 0.3], [0.3, 0.28]])
    @test convex_hull(points) == to_N(N, [[0.1, 0.3], [0.2, 0.1], [0.9, 0.2], [0.4, 0.6]])

    # ===============
    # ConvexHullArray
    # ===============

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(ConvexHull) == ConvexHullArray
    @test LazySets.is_array_constructor(ConvexHullArray)

    # convex hull array of 2 sets
    cha = ConvexHullArray([b1, b2])
    # constructor with size hint and type
    ConvexHullArray(10, N)
    # test alias
    @test CHArray([b1, b2]) isa ConvexHullArray
    # test dimension
    @test dim(cha) == 2
    # test support vector
    d = N[1, 0]
    @test σ(d, cha) == N[2, 2]
    d = N[-1, 0]
    @test σ(d, cha) == N[-1, 0]
    d = N[0, 1]
    @test σ(d, cha) == N[1, 3]
    d = N[0, -1]
    @test σ(d, cha) == N[0, -1]

    # boundedness
    @test isbounded(cha)
    @test !isbounded(ConvexHullArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))

    # isempty
    @test !isempty(cha)

    # test convex hull array of singleton
    ConvexHullArray(
        [Singleton(to_N(N, [10, 0.5])), Singleton(to_N(N, [1.1, 0.2])),
         Singleton(to_N(N, [1.4, 0.3])), Singleton(to_N(N, [1.7, 0.5])),
         Singleton(to_N(N, [1.4, 0.8]))])

    # array getter
    v = Vector{LazySet{N}}()
    @test array(ConvexHullArray(v)) ≡ v

    # in-place modification
    cha = ConvexHullArray(LazySet{N}[])
    @test ConvexHull!(b1, b1) isa ConvexHull && length(array(cha)) == 0
    res = ConvexHull!(b1, cha)
    @test res isa ConvexHullArray && length(array(cha)) == 1
    res = ConvexHull!(cha, b1)
    @test res isa ConvexHullArray && length(array(cha)) == 2
    res = ConvexHull!(cha, cha)
    @test res isa ConvexHullArray && length(array(cha)) == 4

    # concatenation of two convex hull arrays
    @test CH(cha, cha) isa ConvexHull

    # ================
    # common functions
    # ================

    # neutral element
    e = EmptySet{N}()
    @test neutral(ConvexHull) == neutral(ConvexHullArray) == EmptySet
    @test CH(b1, e) == CH(e, b1) == b1
    @test CH(cha, e) == CH(e, cha) == cha
end
