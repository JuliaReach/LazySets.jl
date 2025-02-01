for N in [Float64, Rational{Int}, Float32]
    # ConvexHull of two 2D Ball1
    b1 = Ball1(N[0, 0], N(1))
    b2 = Ball1(N[1, 2], N(1))
    # Test Construction
    ch = ConvexHull(b1, b2)
    @test ch == CH(b1, b2)

    # Universe is absorbing
    U = Universe{N}(2)
    @test CH(b1, U) == CH(U, b1) == CH(U, U) == U
    cha = ConvexHullArray([b1, N(2) * b1, N(3) * b1])
    @test CH(cha, U) == CH(U, cha) == U

    # swap
    ch2 = swap(ch)
    @test ch.X == ch2.Y && ch.Y == ch2.X

    # flatten
    b3 = BallInf(N[0, 0], N(1))
    for M3 in (ConvexHull(ConvexHull(b1, ConvexHullArray([b2])), b3),
               ConvexHullArray([ConvexHull(b1, ConvexHullArray([b2])), b3]))
        M3f = flatten(M3)
        @test M3f isa ConvexHullArray && array(M3f) == [b1, b2, b3]
    end

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
    @test isbounded(ch) && isboundedtype(typeof(ch))
    ch2 = CH(Singleton(N[1]), HalfSpace(N[1], N(1)))
    @test !isbounded(ch2) && !isboundedtype(typeof(ch2))

    # ispolyhedral
    @test ispolyhedral(ch)
    if N isa AbstractFloat
        ch2 = CH(b1, Ball2(N[0, 0], N(1)))
        @test !ispolyhedral(ch2)
    end

    # isempty
    @test !isempty(ch)

    # concretize
    @test LazySets.concrete_function(ConvexHull) == convex_hull
    @test concretize(ch) == convex_hull(ch.X, ch.Y)

    # ===============
    # ConvexHullArray
    # ===============

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(ConvexHull) == ConvexHullArray
    @test LazySets.binary_constructor(ConvexHullArray) == ConvexHull
    @test !LazySets.is_array_constructor(ConvexHull)
    @test LazySets.is_array_constructor(ConvexHullArray)

    # convex hull array of 2 sets
    cha = ConvexHullArray([b1, b2])

    # array interface
    @test array(ch) == array(cha) == [b1, b2]
    @test ch[1] == cha[1] == b1
    @test ch[1:2] == cha[1:2] == [b1, b2]
    @test ch[end] == cha[end] == b2
    @test length(ch) == length(cha) == 2
    v = Vector{LazySet{N}}()
    @test array(ConvexHullArray(v)) ≡ v

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

    # getindex & length
    @test cha[1] == b1 && cha[2] == b2
    @test length(cha) == 2

    # boundedness
    @test isbounded(cha) && isboundedtype(typeof(cha))
    cha2 = ConvexHullArray([Singleton(N[1]), HalfSpace(N[1], N(1))])
    @test !isbounded(cha2) && !isboundedtype(typeof(cha2))

    # isempty
    @test !isempty(cha)

    # test convex hull array of singleton
    ConvexHullArray([Singleton(N[10, 1 // 2]), Singleton(N[11 // 10, 1 // 5]),
                     Singleton(N[7 // 5, 3 // 10]), Singleton(N[17 // 10, 1 // 2]),
                     Singleton(N[7 // 5, 4 // 5])])

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

    # list of vertices of the convex hull array of singleton
    a = Singleton(N[-1, 0])
    b = Singleton(N[1, 0])
    c = Singleton(N[0, 1])
    Δ = ConvexHullArray([a, b, c])
    vlist = vertices_list(Δ)
    @test ispermutation(vlist, [N[-1, 0], N[1, 0], N[0, 1]])

    vΔ = convert(VPolytope, Δ)
    vlist = vertices_list(vΔ)
    @test ispermutation(vlist, [N[-1, 0], N[1, 0], N[0, 1]])

    Δbox = overapproximate(Δ, Hyperrectangle)
    @test Δ ⊆ Δbox && !(Δbox ⊆ Δ)

    # concretize
    cha2 = ConvexHullArray([b1, b2])
    @test concretize(cha2) == convex_hull(b1, b2)

    # ================
    # common functions
    # ================

    # neutral element
    e = EmptySet{N}(2)
    @test neutral(ConvexHull) == neutral(ConvexHullArray) == EmptySet
    @test CH(b1, e) == CH(e, b1) == b1
    @test CH(cha, e) == CH(e, cha) == cha
    @test CH(e, e) == e

    # vertices list
    B1 = BallInf(zeros(N, 2), N(1))
    B2 = Ball1(ones(N, 2), N(1))
    Chull = ConvexHull(B1, B2)
    CHarr = ConvexHullArray([B1, B2])
    for X in [Chull, CHarr]
        @test ispermutation(vertices_list(X),
                            [N[1, -1], N[-1, 1], N[-1, -1], N[1, 2], N[2, 1]])
        @test ispermutation(vertices_list(X; apply_convex_hull=false),
                            [N[1, 1], N[1, -1], N[-1, 1], N[-1, -1],
                             N[1, 2], N[1, 0], N[2, 1], N[0, 1]])
    end

    # linear_combination with empty argument
    for X in (Chull, CHarr)
        for Y in (linear_combination(X, e), linear_combination(e, X))
            @test isempty(Y) && Y isa LazySet{N} && dim(Y) == 2
        end
    end
end
