for N in @tN([Float64, Float32, Rational{Int}])
    # Sum of 2D centered balls in norm 1 and infinity
    b1 = BallInf(N[0, 0], N(2))
    b2 = Ball1(N[0, 0], N(1))
    # Test Construction
    ms = MinkowskiSum(b1, b2)
    @test ms.X == b1
    @test ms.Y == b2

    # convenience constructors
    @test b1 + b2 == b1 ⊕ b2 == ms
    msa = MinkowskiSumArray([b1, b2, b1])
    @test b1 + b2 + b1 == +(b1, b2, b1) == ⊕(b1, b2, b1) == +([b1, b2, b1]) == ⊕([b1, b2, b1]) ==
          msa
    @test +(b1) == ⊕(b1) == b1

    # Universe is absorbing
    # U = Universe{N}(2)  # TODO requires #1099
    # @test b1 ⊕ U == U ⊕ b1 == U ⊕ U == U
    # msa = MinkowskiSumArray([b1, N(2) * b1, N(3) * b1])
    # @test msa ⊕ U == U ⊕ msa == U

    # array interface
    @test array(ms) == [b1, b2] && array(msa) == [b1, b2, b1]
    @test ms[1] == msa[1] == b1
    @test ms[1:2] == msa[1:2] == [b1, b2]
    @test ms[end] == b2 && msa[end] == b1
    @test length(ms) == 2 && length(msa) == 3
    v = Vector{LazySet{N}}()
    @test array(MinkowskiSumArray(v)) ≡ v

    # swap
    ms2 = swap(ms)
    @test ms.X == ms2.Y && ms.Y == ms2.X

    # flatten
    b3 = Ball1(N[1, 1], N(1))
    for M3 in (MinkowskiSum(MinkowskiSum(b1, MinkowskiSumArray([b2])), b3),
               MinkowskiSumArray([MinkowskiSum(b1, MinkowskiSumArray([b2])), b3]))
        M3f = flatten(M3)
        @test M3f isa MinkowskiSumArray && array(M3f) == [b1, b2, b3]
    end

    # Test Dimension
    @test dim(ms) == 2
    # Test Support Vector
    d = N[1, 0]
    v = σ(d, ms)
    @test v[1] == N(3)
    d = N[-1, 0]
    v = σ(d, ms)
    @test v[1] == N(-3)
    d = N[0, 1]
    v = σ(d, ms)
    @test v[2] == N(3)
    d = N[0, -1]
    v = σ(d, ms)
    @test v[2] == N(-3)

    # Sum of not-centered 2D balls in norm 1 and infinity
    b1 = BallInf(N[-1, 3], N(2))
    b2 = Ball1(N[1, 2], N(1))
    s = b1 + b2
    # Test Support Vector
    d = N[1, 0]
    v = σ(d, s)
    @test v[1] == N(3)
    d = N[-1, 0]
    v = σ(d, s)
    @test v[1] == N(-3)
    d = N[0, 1]
    v = σ(d, s)
    @test v[2] == N(8)
    d = N[0, -1]
    v = σ(d, s)
    @test v[2] == N(2)

    # center
    @test center(s) == N[0, 5]

    # support function
    # 2-elements
    ms = MinkowskiSum(Singleton(N[1]), Singleton(N[2]))
    @test ρ(N[1], ms) == 3
    @test ρ(N[-1], ms) == -3
    # 3-elements
    ms = MinkowskiSum(Singleton(N[1]), MinkowskiSum(Singleton(N[2]), Singleton(N[3])))
    @test ρ(N[1], ms) == N(6)
    @test ρ(N[-1], ms) == N(-6)

    # boundedness
    @test isbounded(ms) && isboundedtype(typeof(ms))
    ms2 = Singleton(N[1]) + HalfSpace(N[1], N(1))
    @test !isbounded(ms2) && !isboundedtype(typeof(ms2))

    # ispolyhedral
    @test ispolyhedral(ms)
    if N isa AbstractFloat
        ms2 = MinkowskiSum(b1, Ball2(N[0, 0], N(1)))
        @test !ispolyhedral(ms2)
    end

    # isempty
    @test !isempty(ms)

    # an_element function (falls back to the default implementation
    X = MinkowskiSum(LineSegment(N[0, 0], N[1, 1]),
                     LineSegment(N[1, 0], N[0, 1]))
    v = an_element(X)
    @test v ∈ Ball1(N[1, 1], N(1))

    # constraints list
    B = BallInf(N[0, 0], N(1))
    ms = MinkowskiSum(B, B)
    clist = constraints_list(ms)
    P = HPolytope(clist)
    H = Hyperrectangle(N[0, 0], N[2, 2])
    @test length(clist) == 4

    # membership in the sum of a singleton and a polytopic set
    p = N[1, 2]
    S = Singleton(N[1, 2])
    X = BallInf(zeros(N, 2), N(4))
    @test p ∈ S ⊕ X

    # vertices for minkowski sum of zonotopic sets
    X = Zonotope(N[0, 0], [N[1, 0]])
    Y = Zonotope(N[0, 0], [N[1, 1]])
    @test ispermutation(vertices_list(X + Y), [N[2, 1], N[0, 1], N[-2, -1], N[0, -1]])

    # concretize
    @test LazySets.concrete_function(MinkowskiSum) == minkowski_sum
    @test concretize(ms) == minkowski_sum(B, B)

    # linear_map
    M = N[1 0; 0 1]
    Z = linear_map(M, X + Y)  # identity
    @test isequivalent(X + Y, Z)

    # tohrep/tovrep
    X = Singleton(N[1, 1])
    Y = Singleton(N[2, 3])
    Z = tohrep(X + Y)
    @test Z isa HPolyhedron{N} &&
          ispermutation(constraints_list(Z), constraints_list(X + Y))
    @test tovrep(X + Y) == VPolytope([N[3, 4]])

    # =================
    # MinkowskiSumArray
    # =================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(MinkowskiSum) == MinkowskiSumArray
    @test LazySets.binary_constructor(MinkowskiSumArray) == MinkowskiSum
    @test !LazySets.is_array_constructor(MinkowskiSum)
    @test LazySets.is_array_constructor(MinkowskiSumArray)

    # constructor with size hint and type
    MinkowskiSumArray(10, N)

    # in-place modification
    msa = MinkowskiSumArray(LazySet{N}[])
    @test MinkowskiSum!(b1, b1) isa MinkowskiSum && length(array(msa)) == 0
    res = MinkowskiSum!(b1, msa)
    @test res isa MinkowskiSumArray && length(array(msa)) == 1
    res = MinkowskiSum!(msa, b1)
    @test res isa MinkowskiSumArray && length(array(msa)) == 2
    res = MinkowskiSum!(msa, msa)
    @test res isa MinkowskiSumArray && length(array(msa)) == 4

    ms = MinkowskiSum(b1, b2)
    msa = MinkowskiSumArray([b1, b2])

    # getindex & length
    @test msa[1] == b1 && msa[2] == b2
    @test length(msa) == 2

    # dimension
    @test dim(msa) == 2

    # support vector
    d = N[1, 1]
    @test σ(d, ms) == σ(d, msa)
    d = N[-1, 1]
    @test σ(d, ms) == σ(d, msa)

    # boundedness
    @test isbounded(msa) && isboundedtype(typeof(msa))
    msa2 = MinkowskiSumArray([Singleton(N[1]), HalfSpace(N[1], N(1))])
    @test !isbounded(msa2) && !isboundedtype(typeof(msa2))

    # isempty
    @test !isempty(msa)

    # center
    @test center(msa) == N[0, 5]

    # convert m-sum of m-sum array to m-sum array (#1678)
    s = Singleton(N[1, 2, 3])
    X = s ⊕ MinkowskiSumArray([s, s])
    @test convert(MinkowskiSumArray, X) isa MinkowskiSumArray{N,Singleton{N,Vector{N}}}

    # concretize
    msa = MinkowskiSumArray([B, B])
    @test concretize(msa) == minkowski_sum(B, B)

    # =================
    # CachedMinkowskiSumArray
    # =================

    # caching Minkowski sum
    cms = CachedMinkowskiSumArray(2, N)
    x1 = BallInf(ones(N, 3), N(3))
    x2 = Ball1(ones(N, 3), N(5))
    d = ones(N, 3)
    a = array(cms)
    push!(a, x1)
    σ(d, cms)
    push!(a, x2)
    svec = σ(d, cms)
    @test σ(d, cms) == svec
    @test dim(cms) == 3
    idx = forget_sets!(cms)
    @test idx == 2 && isempty(array(cms)) && cms.cache[d][1] == 0
    push!(a, x1)
    σ(d, cms)
    push!(a, x2)
    idx = forget_sets!(cms)
    @test idx == 1 && length(array(cms)) == 1 && cms.cache[d][1] == 0
    # test issue #367: inplace modification of the direction modified the cache
    d[1] = zero(N)
    @test haskey(cms.cache, ones(N, 3))
    # getindex
    cp = LazySets.CachedPair(1, N[2])
    @test cp[1] == 1 && cp[2] == N[2]
    @test_throws ErrorException cp[3]

    # boundedness
    @test isbounded(cms)
    @test !isbounded(CachedMinkowskiSumArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))

    # isempty
    @test !isempty(cms)

    # concretize
    cms = CachedMinkowskiSumArray([B, B])
    @test concretize(cms) == minkowski_sum(B, B)

    # ================
    # common functions
    # ================

    # neutral element
    z = ZeroSet{N}(2)
    msa = MinkowskiSumArray(LazySet{N}[])
    cms = CachedMinkowskiSumArray(LazySet{N}[])
    @test neutral(MinkowskiSum) == neutral(MinkowskiSumArray) ==
          neutral(CachedMinkowskiSumArray) == ZeroSet
    @test b1 + z == z + b1 == b1
    @test msa + z == z + msa == msa
    @test cms + z == z + cms == cms
    # absorbing element
    e = EmptySet{N}(2)
    @test absorbing(MinkowskiSum) == absorbing(MinkowskiSumArray) ==
          absorbing(CachedMinkowskiSumArray) == EmptySet
    @test b1 + e == e + b1 == msa + e == e + msa == cms + e == e + cms ==
          e + e == e
    # mix of neutral and absorbing element
    @test z + e == e + z == e
end

for N in [Float64]
    # inclusion
    ms = MinkowskiSum(BallInf(N[0, 0], N(1)), BallInf(N[0, 0], N(1)))
    clist = constraints_list(ms)
    P = HPolytope(clist)
    H = Hyperrectangle(N[0, 0], N[2, 2])
    @test length(clist) == 4 && P ⊆ H && H ⊆ P

    # concrete cartesian product
    if test_suite_polyhedra
        X = BallInf(N[0, 0], N(1)) ⊕ Singleton(N[1, 2])
        Y = Ball1(N[1, 1], N(1))
        A = cartesian_product(X, Y; algorithm="vrep")
        Ah = tohrep(A; backend=CDDLib.Library())
        B = cartesian_product(X, Y; algorithm="hrep")
        @test Ah ⊆ B && B ⊆ Ah
    end
end
