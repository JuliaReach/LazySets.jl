for N in [Float64, Rational{Int}, Float32]
    # Sum of 2D centered balls in norm 1 and infinity
    b1 = BallInf(N[0, 0], N(2))
    b2 = Ball1(N[0, 0], N(1))
    # Test Construction
    X = MinkowskiSum(b1, b2)
    @test X.X == b1
    @test X.Y == b2

    # swap
    ms2 = swap(X)
    @test X.X == ms2.Y && X.Y == ms2.X

    # Test Dimension
    @test dim(X) == 2
    # Test Support Vector
    d = N[1, 0]
    v = σ(d, X)
    @test v[1] == N(3)
    d = N[-1, 0]
    v = σ(d, X)
    @test v[1] == N(-3)
    d = N[0, 1]
    v = σ(d, X)
    @test v[2] == N(3)
    d = N[0, -1]
    v = σ(d, X)
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

    # Sum of array of LazySet
    # 2-elements
    ms = MinkowskiSum(Singleton(N[1]), Singleton(N[2]))
    @test ρ(N[1], ms) == 3
    @test ρ(N[-1], ms) == -3
    # 3-elements
    ms = MinkowskiSum(Singleton(N[1]), MinkowskiSum(Singleton(N[2]), Singleton(N[3])))
    @test ρ(N[1], ms) == N(6)
    @test ρ(N[-1], ms) == N(-6)

    # boundedness
    @test isbounded(ms)
    @test !isbounded(Singleton(N[1]) + HalfSpace(N[1], N(1)))

    # isempty
    @test !isempty(ms)

    # an_element function (falls back to the default implementation
    X = MinkowskiSum(LineSegment(N[0, 0], N[1, 1]),
                     LineSegment(N[1, 0], N[0, 1]))
    v = an_element(X)
    @test v ∈ Ball1(N[1, 1], N(1))

    # constraints list
    ms = MinkowskiSum(BallInf(N[0, 0], N(1)), BallInf(N[0, 0], N(1)))
    clist = constraints_list(ms)
    P = HPolytope(clist)
    H = Hyperrectangle(N[0, 0], N[2, 2])
    @test length(clist) == 4

    # =================
    # MinkowskiSumArray
    # =================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(MinkowskiSum) == MinkowskiSumArray
    @test LazySets.is_array_constructor(MinkowskiSumArray)

    # array getter
    v = Vector{LazySet{N}}()
    @test array(MinkowskiSumArray(v)) ≡ v

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

    # dimension
    @test dim(msa) == 2

    # support vector
    d = N[1, 1]
    @test σ(d, ms) == σ(d, msa)
    d = N[-1, 1]
    @test σ(d, ms) == σ(d, msa)

    # boundedness
    @test isbounded(msa)
    @test !isbounded(MinkowskiSumArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))

    # isempty
    @test !isempty(msa)

    # =================
    # CacheMinkowskiSum
    # =================

    # caching Minkowski sum
    cms = CacheMinkowskiSum(2, N)
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
    @test !isbounded(CacheMinkowskiSum([Singleton(N[1]), HalfSpace(N[1], N(1))]))

    # isempty
    @test !isempty(cms)

    # ================
    # common functions
    # ================

    # neutral element
    z = ZeroSet{N}(2)
    msa = MinkowskiSumArray(LazySet{N}[])
    cms = CacheMinkowskiSum(LazySet{N}[])
    @test neutral(MinkowskiSum) == neutral(MinkowskiSumArray) ==
          neutral(CacheMinkowskiSum) == ZeroSet
    @test b1 + z == z + b1 == b1
    @test msa + z == z + msa == msa
    @test cms + z == z + cms == cms
    # absorbing element
    e = EmptySet{N}()
    @test absorbing(MinkowskiSum) == absorbing(MinkowskiSumArray) ==
          absorbing(CacheMinkowskiSum) == EmptySet
    @test b1 + e == e + b1 == msa + e == e + msa == cms + e == e + cms ==
          e + e == e
    # mix of neutral and absorbing element
    @test z + e == e + z == e
end

# tests that only work with Float64
for N in [Float64]
    # inclusion
    ms = MinkowskiSum(BallInf(N[0, 0], N(1)), BallInf(N[0, 0], N(1)))
    clist = constraints_list(ms)
    P = HPolytope(clist)
    H = Hyperrectangle(N[0, 0], N[2, 2])
    @test length(clist) == 4 && P ⊆ H && H ⊆ P
end
